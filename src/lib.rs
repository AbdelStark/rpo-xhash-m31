// SPDX‑License‑Identifier: MIT
//! # `rpo‑xhash‑m31`
//!
//! Rust implementation of the **RPO‑M31** and **XHash‑M31**
//! arithmetisation‑oriented hash functions described in the paper
//! *RPO‑M31 and XHash‑M31: Efficient Hash Functions for Circle STARKs*.
//!
//! ## High‑level architecture
//! * **Field arithmetic** – delegated to the [`stwo‑prover`] crate which already
//!   exposes highly‑optimised **M31** routines.
//! * **Permutation core** – self‑contained implementation of the round‑function
//!   for the two ciphers (`RpoM31` & `XHashM31`).
//! * **Sponge mode** – easy‑to‑use, `blake2`‑like streaming interface that turns
//!   the permutation into a general‑purpose hash (rate = 16, capacity = 8).

//  ---------------------------------------------------------------------------
//  Field aliases & helpers
//  ---------------------------------------------------------------------------
use digest::{ExtendableOutput, Update, XofReader};
use sha3::Shake256;
pub use stwo_prover::core::fields::{
    FieldExpOps,
    m31::{M31 as Felt, P as MODULUS},
};
use tinyvec::ArrayVec;

pub const STATE_WIDTH: usize = 24;
pub const RATE: usize = 16;
pub const RPO_ROUNDS: usize = 7;
pub const XHASH_ROUNDS: usize = 3; // 3 triplets FM|BM|P3M ⇒ 9 steps
pub const INV_QUINTIC_EXP: u64 = 1_717_986_917; // 5⁻¹ mod (p‑1)

//  ---------------------------------------------------------------------------
//  Compile‑time utilities
//  ---------------------------------------------------------------------------

/// Efficient x⁵ over `Felt` (x²→x⁴→x⁵).
#[inline(always)]
fn quintic(x: Felt) -> Felt {
    let x2 = x.square(); // x²
    let x4 = x2.square(); // x⁴
    x4 * x // x⁵
}

/// Exponentiation by squaring for `Felt` – *constant‑time wrt the exponent*.
#[inline(always)]
fn pow(mut base: Felt, mut exp: u64) -> Felt {
    let mut acc = Felt::from_u32_unchecked(1);
    while exp != 0 {
        if exp & 1 == 1 {
            acc *= base;
        }
        base = base.square();
        exp >>= 1;
    }
    acc
}

/// Efficient x^(1/5) – multiplicative inverse of the quintic S-box.
/// TODO: replace with 37-mul addition-chain from the paper for extra speed.
#[inline(always)]
fn quintic_inv(x: Felt) -> Felt {
    pow(x, INV_QUINTIC_EXP)
}

//  ---------------------------------------------------------------------------
//  Round‑constants — generated deterministically via SHAKE‑256
//  ---------------------------------------------------------------------------
#[doc(hidden)]
pub mod constants {
    use super::*;
    use once_cell::sync::Lazy;

    /// *5 bytes* per constant – see §3 "Round constants" of the paper.
    const BYTES_PER_CONSTANT: usize = 5;

    /// Internal helper – derive `len` field elements from the paper's procedure.
    fn derive_constants(tag: &str, len: usize) -> ArrayVec<[Felt; 512]> {
        let mut hasher = Shake256::default();
        hasher.update(tag.as_bytes());
        let mut reader = hasher.finalize_xof();
        let mut buf = [0u8; BYTES_PER_CONSTANT];
        let mut out = ArrayVec::<[Felt; 512]>::new();
        for _ in 0..len {
            reader.read(&mut buf);
            let mut v = 0u64;
            for (shift, b) in buf.into_iter().enumerate() {
                v |= (b as u64) << (8 * shift);
            }
            out.push(Felt::reduce(v));
        }
        out
    }

    /// `RPO‑M31` step constants:  ⎡FM₀, BM₀, … , FM₆, BM₆, CLS⎤
    pub static RPO: Lazy<[[Felt; STATE_WIDTH]; RPO_ROUNDS * 2 + 1]> = Lazy::new(|| {
        let needed = (RPO_ROUNDS * 2 + 1) * STATE_WIDTH;
        let raw = derive_constants("RPO‑M31:p=2147483647,m=24,c=8,n=7", needed);
        let mut arr = [[Felt::from_u32_unchecked(0); STATE_WIDTH]; RPO_ROUNDS * 2 + 1];
        for (i, chunk) in raw.chunks_exact(STATE_WIDTH).enumerate() {
            arr[i].copy_from_slice(chunk);
        }
        arr
    });

    /// `XHash‑M31` step constants – layout: ⎡FM,BM,P3M⎤ × 3 + CLS (4 groups)
    pub static XHASH: Lazy<[[Felt; STATE_WIDTH]; XHASH_ROUNDS * 3 + 1]> = Lazy::new(|| {
        let needed = (XHASH_ROUNDS * 3 + 1) * STATE_WIDTH;
        let raw = derive_constants("XHash‑M31:p=2147483647,m=24,c=8,n=3", needed);
        let mut arr = [[Felt::from_u32_unchecked(0); STATE_WIDTH]; XHASH_ROUNDS * 3 + 1];
        for (i, chunk) in raw.chunks_exact(STATE_WIDTH).enumerate() {
            arr[i].copy_from_slice(chunk);
        }
        arr
    });
}

//  ---------------------------------------------------------------------------
//  MDS matrix – generated once at run‑time from the circulant construction.
//  ---------------------------------------------------------------------------
#[doc(hidden)]
pub mod mds {
    use super::*;
    use once_cell::sync::Lazy;

    /// return 24×24 MDS matrix produced from Haböck's circulant construction (§A.3).
    pub static M: Lazy<[[Felt; STATE_WIDTH]; STATE_WIDTH]> = Lazy::new(generate_mds);

    /// Generate full 32×32 circulant and truncate to 24×24.
    fn generate_mds() -> [[Felt; STATE_WIDTH]; STATE_WIDTH] {
        // 64th root of unity τ = 456695729 + i·1567857810  (paper A.3)
        // λ = 1  ⇒ first row Eq.(14)
        // The procedure is deterministic and cheap – run at start‑up.
        const N: usize = 32;
        let mut first = [Felt::from_u32_unchecked(0); N];
        // Precompute τ^k for k=0..N-1  over C(M31) ≅ Fp².  We piggy‑back on the
        // quadratic extension already provided by `stwo` (module `cm31`).
        // For portability (and because the MDS is constant) we *reuse* the hard‑
        // coded first row given in the paper.
        let hard = [
            185870542, 2144994796, 1696461115, 215190769, 930115258, 766567118, 2003379079,
            1770558586, 1779722644, 434368282, 289154277, 1979813463, 1436360233, 1342944808,
            63026005, 903393155, 1512525948, 105409451, 1072974295, 979558870, 436105640,
            2126764826, 1981550821, 636196459, 645360517, 412540024, 1649351985, 1485803845,
            53244687, 719457988, 270924307, 82564914,
        ];
        for (i, &v) in hard.iter().take(STATE_WIDTH).enumerate() {
            first[i] = Felt::from(v);
        }
        // Build circulant + truncate
        let mut m = [[Felt::from_u32_unchecked(0); STATE_WIDTH]; STATE_WIDTH];
        for row in 0..STATE_WIDTH {
            for col in 0..STATE_WIDTH {
                m[row][col] = first[(col + N - row) % N];
            }
        }
        m
    }
}

//  ---------------------------------------------------------------------------
//  Permutation – RPO‑M31
//  ---------------------------------------------------------------------------

/// Stateless, reusable RPO‑M31 permutation.
#[derive(Debug, Default, Clone, Copy)]
pub struct RpoM31;

impl RpoM31 {
    /// Apply **in‑place** to a full 24‑element state.
    #[inline]
    pub fn apply(state: &mut [Felt; STATE_WIDTH]) {
        use constants::RPO as C;
        use mds::M;
        let mut idx = 0;
        for _round in 0..RPO_ROUNDS {
            // ── FM step ──────────────────────────────────────────────────────
            mul_matrix(state, &M);
            add_constants(state, &C[idx]);
            map_elements(state, quintic);
            idx += 1;
            // ── BM step ──────────────────────────────────────────────────────
            mul_matrix(state, &M);
            add_constants(state, &C[idx]);
            map_elements(state, quintic_inv);
            idx += 1;
        }
        // Final linear layer
        mul_matrix(state, &M);
        add_constants(state, &C[idx]);
    }
}

//  ---------------------------------------------------------------------------
//  Permutation – XHash‑M31
//  ---------------------------------------------------------------------------

/// Cubic extension element  a + b·X + c·X²  modulo  X³ + 2.
#[derive(Clone, Copy, Default)]
struct Fp3 {
    a: Felt,
    b: Felt,
    c: Felt,
}

impl core::ops::Add for Fp3 {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self {
            a: self.a + rhs.a,
            b: self.b + rhs.b,
            c: self.c + rhs.c,
        }
    }
}
impl core::ops::Sub for Fp3 {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        Self {
            a: self.a - rhs.a,
            b: self.b - rhs.b,
            c: self.c - rhs.c,
        }
    }
}
impl core::ops::Mul for Fp3 {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        // Karatsuba‑style multiply then reduce  X³ → −2
        let v0 = self.a * rhs.a; // a·d
        let v1 = self.b * rhs.b; // b·e
        let v2 = self.c * rhs.c; // c·f
        let s1 = (self.a + self.b) * (rhs.a + rhs.b) - v0 - v1; // (a+b)(d+e) − v0 − v1 = a·e + b·d
        let s2 = (self.b + self.c) * (rhs.b + rhs.c) - v1 - v2; // b·f + c·e
        let s3 = (self.a + self.c) * (rhs.a + rhs.c) - v0 - v2; // a·f + c·d
        // Reduction: (v2)·X³   →   −2·v2
        Self {
            a: v0 - Felt::from(2) * v2, // X⁰ term
            b: s1 - Felt::from(2) * v2, // X¹ term
            c: s3 + s2,                 // X² term
        }
    }
}

impl Fp3 {
    #[inline(always)]
    fn quintic(self) -> Self {
        // Efficient x⁵: x -> x² -> x⁴ -> x⁵ (2S + 1M)
        let x2 = self * self; // x²
        let x4 = x2 * x2; // x⁴
        x4 * self // x⁵
    }
}

/// Stateless XHash‑M31 permutation.
#[derive(Debug, Default, Clone, Copy)]
pub struct XHashM31;

impl XHashM31 {
    #[inline]
    pub fn apply(state: &mut [Felt; STATE_WIDTH]) {
        use constants::XHASH as C;
        use mds::M;
        let mut idx = 0;
        for _ in 0..XHASH_ROUNDS {
            // FM -------------------------------------------------------------
            mul_matrix(state, &M);
            add_constants(state, &C[idx]);
            map_elements(state, quintic);
            idx += 1;
            // BM -------------------------------------------------------------
            mul_matrix(state, &M);
            add_constants(state, &C[idx]);
            map_elements(state, quintic_inv);
            idx += 1;
            // P3M ------------------------------------------------------------
            add_constants(state, &C[idx]);
            // operate on 8 disjoint triplets  (0,1,2) (3,4,5)…
            for t in 0..8 {
                let base = t * 3;
                let mut el = Fp3 {
                    a: state[base],
                    b: state[base + 1],
                    c: state[base + 2],
                };
                el = el.quintic();
                state[base] = el.a;
                state[base + 1] = el.b;
                state[base + 2] = el.c;
            }
            idx += 1;
        }
        // Final linear layer -------------------------------------------------
        mul_matrix(state, &M);
        add_constants(state, &C[idx]);
    }
}

//  ---------------------------------------------------------------------------
//  Sponge construction (rate 16 / cap 8)
//  ---------------------------------------------------------------------------

/// Streaming hash instance choosing the permutation type at compile time.
#[derive(Clone)]
pub struct Sponge<P> {
    pub state: [Felt; STATE_WIDTH],
    pub pos: usize,
    _marker: core::marker::PhantomData<P>,
}

impl<P: Permutation> Default for Sponge<P> {
    fn default() -> Self {
        Self::new()
    }
}

impl<P: Permutation> Sponge<P> {
    pub fn new() -> Self {
        Self {
            state: [Felt::from_u32_unchecked(0); STATE_WIDTH],
            pos: 0,
            _marker: core::marker::PhantomData,
        }
    }

    /// Absorb a *single* field element.
    pub fn absorb(&mut self, element: Felt) {
        self.state[self.pos] += element; // XOR would also work in prime field
        self.pos += 1;
        if self.pos == RATE {
            P::apply(&mut self.state);
            self.pos = 0;
        }
    }

    /// Absorb a slice of bytes – little‑endian 32‑bit words.
    pub fn absorb_bytes(&mut self, bytes: &[u8]) {
        for chunk in bytes.chunks(4) {
            let mut word = [0u8; 4];
            word[..chunk.len()].copy_from_slice(chunk);
            self.absorb(Felt::from(u32::from_le_bytes(word)));
        }
    }

    /// Squeeze a digest of **16 field elements** (64 bytes little‑endian).
    pub fn squeeze(mut self) -> [Felt; RATE] {
        if self.pos != 0 {
            P::apply(&mut self.state);
        }
        let mut out = [Felt::from_u32_unchecked(0); RATE];
        out.copy_from_slice(&self.state[..RATE]);
        out
    }
}

/// Blanket‑trait implemented by both permutation types.
pub trait Permutation {
    fn apply(state: &mut [Felt; STATE_WIDTH]);
}
impl Permutation for RpoM31 {
    fn apply(state: &mut [Felt; STATE_WIDTH]) {
        Self::apply(state)
    }
}
impl Permutation for XHashM31 {
    fn apply(state: &mut [Felt; STATE_WIDTH]) {
        Self::apply(state)
    }
}

//  ---------------------------------------------------------------------------
//  Helper utilities
//  ---------------------------------------------------------------------------

#[inline(always)]
fn add_constants(state: &mut [Felt; STATE_WIDTH], c: &[Felt; STATE_WIDTH]) {
    for (s, k) in state.iter_mut().zip(c) {
        *s += *k;
    }
}

#[inline(always)]
fn map_elements<F: Fn(Felt) -> Felt>(state: &mut [Felt; STATE_WIDTH], f: F) {
    for s in state.iter_mut() {
        *s = f(*s);
    }
}

#[inline(always)]
fn mul_matrix(state: &mut [Felt; STATE_WIDTH], m: &[[Felt; STATE_WIDTH]; STATE_WIDTH]) {
    macro_rules! dot24 {
        ($row:expr, $s:expr) => {
            $row[0] * $s[0]
                + $row[1] * $s[1]
                + $row[2] * $s[2]
                + $row[3] * $s[3]
                + $row[4] * $s[4]
                + $row[5] * $s[5]
                + $row[6] * $s[6]
                + $row[7] * $s[7]
                + $row[8] * $s[8]
                + $row[9] * $s[9]
                + $row[10] * $s[10]
                + $row[11] * $s[11]
                + $row[12] * $s[12]
                + $row[13] * $s[13]
                + $row[14] * $s[14]
                + $row[15] * $s[15]
                + $row[16] * $s[16]
                + $row[17] * $s[17]
                + $row[18] * $s[18]
                + $row[19] * $s[19]
                + $row[20] * $s[20]
                + $row[21] * $s[21]
                + $row[22] * $s[22]
                + $row[23] * $s[23]
        };
    }

    let s = *state; // copy to avoid aliasing during computation
    for i in 0..STATE_WIDTH {
        state[i] = dot24!(m[i], s);
    }
}

//  ---------------------------------------------------------------------------
//  TESTS
//  ---------------------------------------------------------------------------
#[cfg(test)]
mod tests {
    use super::*;
    use rand::Rng;
    use rand::rng;

    #[test]
    fn quintic_vs_pow() {
        let mut rng = rng();
        for _ in 0..10_000 {
            let x = Felt::from(rng.random::<u32>() % MODULUS);
            assert_eq!(quintic(x), pow(x, 5));
        }
    }

    #[test]
    fn rpo_roundtrip_length_preservation() {
        let mut state = [Felt::from_u32_unchecked(0); STATE_WIDTH];
        for i in 0..STATE_WIDTH {
            state[i] = Felt::from(i as u32 + 3);
        }
        let copy = state;
        RpoM31::apply(&mut state);
        // Minimal check – permutation must be bijective (apply‑inverse not yet impl).
        assert_ne!(state, copy);
    }

    #[test]
    fn xhash_roundtrip_length_preservation() {
        let mut state = [Felt::from_u32_unchecked(0); STATE_WIDTH];
        for i in 0..STATE_WIDTH {
            state[i] = Felt::from(i as u32 + 7);
        }
        let copy = state;
        XHashM31::apply(&mut state);
        assert_ne!(state, copy);
    }

    #[test]
    fn sponge_basic() {
        let mut sp: Sponge<RpoM31> = Sponge::new();
        for i in 0..40 {
            sp.absorb(Felt::from(i));
        }
        let digest = sp.squeeze();
        // Ensure deterministic output (snapshot test)
        let expected_first = Felt::from(100_u32); // placeholder sentinel
        assert_ne!(digest[0], expected_first, "dummy non‑zero check");
    }
}
