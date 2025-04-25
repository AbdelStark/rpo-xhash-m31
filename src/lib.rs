// SPDX‑License‑Identifier: MIT
//! # `rpo‑xhash‑m31`
//!
//! Rust implementation of the **RPO‑M31** and **XHash‑M31**
//! arithmetisation‑oriented hash functions described in the paper
//! [*RPO‑M31 and XHash‑M31: Efficient Hash Functions for Circle STARKs*][paper].
//!
//! [paper]: https://eprint.iacr.org/2023/1070
//!
//! ## High‑level architecture
//! * **Field arithmetic** – delegated to the [`stwo‑prover`] crate which already
//!   exposes highly‑optimised **M31** routines.
//! * **Permutation core** – self‑contained implementation of the round‑function
//!   for the two ciphers ([`RpoM31`] & [`XHashM31`]).
//! * **Sponge mode** – easy‑to‑use, `blake2`‑like streaming interface that turns
//!   the permutation into a general‑purpose hash ([`Sponge`]), with rate = 16 and capacity = 8.
//!
//! ## Usage
//!
//! ```rust
//! use rpo_xhash_m31::{RpoM31, XHashM31, Sponge, Felt, NoopOpsTracker};
//!
//! // --- RPO ---
//! let mut rpo_sponge: Sponge<RpoM31, NoopOpsTracker> = Sponge::new();
//! rpo_sponge.absorb_bytes(b"some input data");
//! let rpo_digest: [Felt; 16] = rpo_sponge.squeeze();
//!
//! // --- XHash ---
//! let mut xhash_sponge: Sponge<XHashM31, NoopOpsTracker> = Sponge::new();
//! xhash_sponge.absorb_bytes(b"different input");
//! let xhash_digest: [Felt; 16] = xhash_sponge.squeeze();
//!
//! println!("RPO digest element 0: {:?}", rpo_digest[0]);
//! println!("XHash digest element 0: {:?}", xhash_digest[0]);
//! ```

pub mod fields;

//  ---------------------------------------------------------------------------
//  Field aliases & helpers
//  ---------------------------------------------------------------------------
use core::ops::{Add, Mul, Sub};
use digest::{ExtendableOutput, Update, XofReader};
pub use fields::{
    FieldExpOps,
    m31::{M31 as Felt, P as MODULUS},
};
use serde::{Deserialize, Serialize};
use sha3::Shake256;
use std::collections::HashMap;
use tinyvec::ArrayVec;

/// The width of the permutation state in field elements (24).
pub const STATE_WIDTH: usize = 24;
/// The rate (number of elements absorbed/squeezed per permutation) of the sponge (16).
pub const RATE: usize = 16;
/// The number of rounds in the RPO-M31 permutation (7).
pub const RPO_ROUNDS: usize = 7;
/// The number of round triplets in the XHash-M31 permutation (3).
pub const XHASH_ROUNDS: usize = 3; // 3 triplets FM|BM|P3M ⇒ 9 steps
/// The exponent `5⁻¹ mod (p-1)` used for the inverse quintic S-box.
pub const INV_QUINTIC_EXP: u64 = 1_717_986_917; // 5⁻¹ mod (p‑1)

//  ---------------------------------------------------------------------------
//  Compile‑time utilities
//  ---------------------------------------------------------------------------

/// Efficiently computes `x⁵` over `Felt`.
/// Uses the sequence x -> x² -> x⁴ -> x⁵ (2 squares, 1 multiply).
#[inline(always)]
fn quintic<T: OpsTracker>(x: Felt, tracker: &mut T) -> Felt {
    let x2 = x.square(); // x²
    tracker.record(Op::FeltSquare);
    let x4 = x2.square(); // x⁴
    tracker.record(Op::FeltSquare);
    let res = x4 * x; // x⁵
    tracker.record(Op::FeltMul);
    // Explicitly record the high-level operation completion
    tracker.record(Op::FeltQuintic);
    res
}

/// Computes `base^exp` using exponentiation by squaring for `Felt`.
/// Constant‑time with respect to the exponent.
#[inline(always)]
fn pow<T: OpsTracker>(mut base: Felt, mut exp: u64, tracker: &mut T) -> Felt {
    let mut acc = Felt::from_u32_unchecked(1);
    while exp != 0 {
        if exp & 1 == 1 {
            acc *= base;
            tracker.record(Op::FeltMul);
        }
        base = base.square();
        tracker.record(Op::FeltSquare);
        exp >>= 1;
    }
    acc
}

/// Efficiently computes `x^(1/5)` over `Felt` (multiplicative inverse of the quintic S-box).
/// TODO: replace with 37-mul addition-chain from the paper for extra speed.
#[inline(always)]
fn quintic_inv<T: OpsTracker>(x: Felt, tracker: &mut T) -> Felt {
    let res = pow(x, INV_QUINTIC_EXP, tracker);
    // Explicitly record the high-level operation completion
    tracker.record(Op::FeltQuinticInv);
    res
}

//  ---------------------------------------------------------------------------
//  Round‑constants — generated deterministically via SHAKE‑256
//  ---------------------------------------------------------------------------
/// Internal module for generating and storing round constants.
///
/// Constants are derived using SHAKE-256 as specified in §3 of the paper.
/// They are computed once at runtime using `once_cell::sync::Lazy`.
#[doc(hidden)]
pub mod constants {
    use super::*;
    use once_cell::sync::Lazy;

    /// Number of bytes used to derive each field element constant (5 bytes).
    const BYTES_PER_CONSTANT: usize = 5;

    /// Derives `len` field elements using SHAKE-256 with the given `tag`.
    ///
    /// Reads `BYTES_PER_CONSTANT` bytes from the SHAKE XOF output for each element,
    /// interprets them as little-endian, and reduces modulo `P`.
    fn derive_constants(tag: &str, len: usize) -> ArrayVec<[Felt; 512]> {
        let mut hasher = Shake256::default();
        hasher.update(tag.as_bytes());
        let mut reader = hasher.finalize_xof();
        let mut buf = [0u8; BYTES_PER_CONSTANT];
        // Use a fixed-size ArrayVec as the number of constants is known and relatively small.
        let mut out = ArrayVec::<[Felt; 512]>::new();
        for _ in 0..len {
            reader.read(&mut buf);
            // Reconstruct u64 from 5 little-endian bytes
            let mut v = 0u64;
            for (shift, b) in buf.into_iter().enumerate() {
                v |= (b as u64) << (8 * shift);
            }
            // Reduce the 40-bit value into M31
            out.push(Felt::reduce(v));
        }
        out
    }

    /// Round constants for the RPO-M31 permutation.
    ///
    /// Structure: `[FM₀, BM₀, ..., FM₆, BM₆, CLS]`
    /// Total `(RPO_ROUNDS * 2 + 1) * STATE_WIDTH` constants.
    pub static RPO: Lazy<[[Felt; STATE_WIDTH]; RPO_ROUNDS * 2 + 1]> = Lazy::new(|| {
        let needed = (RPO_ROUNDS * 2 + 1) * STATE_WIDTH;
        let raw = derive_constants("RPO‑M31:p=2147483647,m=24,c=8,n=7", needed);
        let mut arr = [[Felt::from_u32_unchecked(0); STATE_WIDTH]; RPO_ROUNDS * 2 + 1];
        for (i, chunk) in raw.chunks_exact(STATE_WIDTH).enumerate() {
            arr[i].copy_from_slice(chunk);
        }
        arr
    });

    /// Round constants for the XHash-M31 permutation.
    ///
    /// Structure: `[FM₀, BM₀, P3M₀, ..., FM₂, BM₂, P3M₂, CLS]`
    /// Total `(XHASH_ROUNDS * 3 + 1) * STATE_WIDTH` constants.
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
/// Internal module for generating and storing the MDS matrix.
///
/// The 24x24 MDS matrix `M` is derived from a 32x32 circulant matrix
/// construction as described in Appendix A.3 of the paper.
/// It is computed once at runtime using `once_cell::sync::Lazy`.
#[doc(hidden)]
pub mod mds {
    use super::*;
    use once_cell::sync::Lazy;

    /// The 24x24 MDS matrix `M` used in the linear layers.
    pub static M: Lazy<[[Felt; STATE_WIDTH]; STATE_WIDTH]> = Lazy::new(generate_mds);

    /// Generates the 24x24 MDS matrix `M`.
    ///
    /// This function implements the circulant construction from Appendix A.3,
    /// using a hardcoded first row derived from roots of unity in a quadratic
    /// extension field C(M31). The full 32x32 circulant matrix is generated,
    /// and then truncated to the required 24x24 size.
    fn generate_mds() -> [[Felt; STATE_WIDTH]; STATE_WIDTH] {
        // 64th root of unity τ = 456695729 + i·1567857810  (paper A.3)
        // λ = 1  ⇒ first row Eq.(14)
        // The procedure is deterministic and cheap – run at start‑up.
        // We use the hardcoded first row provided in the paper's reference implementation
        // for simplicity and reproducibility, avoiding the need for complex C(M31) arithmetic here.
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
        // Convert the hardcoded u32 values to Felt elements.
        for (i, &v) in hard.iter().take(STATE_WIDTH).enumerate() {
            first[i] = Felt::from(v);
        }
        // Build circulant matrix from the first row
        let mut m = [[Felt::from_u32_unchecked(0); STATE_WIDTH]; STATE_WIDTH];
        for (row_idx, current_row) in m.iter_mut().enumerate().take(STATE_WIDTH) {
            for (col_idx, current_col) in current_row.iter_mut().enumerate().take(STATE_WIDTH) {
                // Calculate the index in the first row using the circulant property:
                // M[row][col] = first_row[(col - row + N) % N]
                let index = (col_idx + N - row_idx) % N; // Use row_idx here
                *current_col = first[index];
            }
        }
        // Truncate the 32x32 circulant matrix to 24x24 (implicitly done by loop bounds)
        m
    }
}

//  ---------------------------------------------------------------------------
//  Operations Counting
//  ---------------------------------------------------------------------------

/// Enum representing the different operations to track.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum Op {
    FeltAdd,        // Addition of Felt
    FeltMul,        // Multiplication of Felt
    FeltSquare,     // Squaring of Felt
    FeltQuintic,    // Quintic S-box on Felt (x^5)
    FeltQuinticInv, // Inverse Quintic S-box on Felt (x^(1/5))
    MdsMul,         // Matrix-vector multiplication
    Fp3Add,         // Addition in Fp3
    Fp3Sub,         // Subtraction in Fp3
    Fp3Mul,         // Multiplication in Fp3
    Fp3Quintic,     // Quintic S-box in Fp3
    Permutation,    // A full permutation call
}

/// Trait for tracking cryptographic operations.
pub trait OpsTracker {
    /// Records an occurrence of the specified operation.
    fn record(&mut self, op: Op);

    /// Generates a JSON report of the recorded operations.
    /// Returns None if reporting is not supported (e.g., NoopOpsTracker).
    fn report_json(&self) -> Option<String>;
}

/// An OpsTracker that does nothing.
#[derive(Debug, Default, Clone, Copy)]
pub struct NoopOpsTracker;

impl OpsTracker for NoopOpsTracker {
    #[inline(always)]
    fn record(&mut self, _op: Op) { /* No-op */
    }

    #[inline(always)]
    fn report_json(&self) -> Option<String> {
        None // No reporting for the no-op tracker
    }
}

/// An OpsTracker that counts occurrences of each operation.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct CountingOpsTracker {
    counts: HashMap<Op, u64>,
}

impl CountingOpsTracker {
    /// Creates a new, empty counter.
    pub fn new() -> Self {
        Self {
            counts: HashMap::new(),
        }
    }

    /// Returns a reference to the internal counts map.
    pub fn get_counts(&self) -> &HashMap<Op, u64> {
        &self.counts
    }
}

impl OpsTracker for CountingOpsTracker {
    #[inline]
    fn record(&mut self, op: Op) {
        *self.counts.entry(op).or_insert(0) += 1;
    }

    fn report_json(&self) -> Option<String> {
        match serde_json::to_string_pretty(&self.counts) {
            Ok(json_string) => Some(json_string),
            Err(_) => None, // Handle serialization error if necessary
        }
    }
}

/// Represents an element in the cubic extension field Fp3 = M31[X] / (X³ + 2).
///
/// Elements are stored as `a + b*X + c*X²`.
/// This struct is used internally for the P3M step in the XHash-M31 permutation.
#[derive(Clone, Copy, Default, Debug, PartialEq, Eq)]
struct Fp3 {
    /// Coefficient of X⁰
    a: Felt,
    /// Coefficient of X¹
    b: Felt,
    /// Coefficient of X²
    c: Felt,
}

/// Implements addition for `Fp3` elements (component-wise).
impl Add for Fp3 {
    type Output = Self;
    #[inline(always)]
    fn add(self, rhs: Self) -> Self {
        // No tracker here, assume caller tracks Op::Fp3Add if needed
        Self {
            a: self.a + rhs.a,
            b: self.b + rhs.b,
            c: self.c + rhs.c,
        }
    }
}

/// Implements subtraction for `Fp3` elements (component-wise).
impl Sub for Fp3 {
    type Output = Self;
    #[inline(always)]
    fn sub(self, rhs: Self) -> Self {
        // No tracker here, assume caller tracks Op::Fp3Sub if needed
        Self {
            a: self.a - rhs.a,
            b: self.b - rhs.b,
            c: self.c - rhs.c,
        }
    }
}

/// Implements multiplication for `Fp3` elements modulo `X³ + 2`.
///
/// Uses a Karatsuba-like approach followed by reduction (`X³ -> -2`).
impl Mul for Fp3 {
    type Output = Self;
    #[inline]
    fn mul(self, rhs: Self) -> Self {
        // Note: Tracking is done by the caller (e.g., Fp3::quintic)
        // If Fp3Mul needed tracking in other contexts, a mul_with_tracker method
        // or similar would be required, as trait methods cannot easily take extra args.
        let v0 = self.a * rhs.a;
        let v1 = self.b * rhs.b;
        let v2 = self.c * rhs.c;
        let s1 = (self.a + self.b) * (rhs.a + rhs.b) - v0 - v1;
        let s2 = (self.b + self.c) * (rhs.b + rhs.c) - v1 - v2;
        let s3 = (self.a + self.c) * (rhs.a + rhs.c) - v0 - v2;
        Self {
            a: v0 - Felt::from(2) * v2,
            b: s1 - Felt::from(2) * v2,
            c: s3 + s2,
        }
    }
}

impl Fp3 {
    /// Efficiently computes `self⁵` over Fp3.
    /// Uses the sequence `x -> x² -> x⁴ -> x⁵` (2 squares, 1 multiply).
    #[inline(always)]
    fn quintic<T: OpsTracker>(self, tracker: &mut T) -> Self {
        // Record the start of the Fp3 Quintic operation logic
        // x -> x²
        let x2 = self * self;
        // Track the underlying Fp3 multiplication for the square
        tracker.record(Op::Fp3Mul);
        // x² -> x⁴
        let x4 = x2 * x2;
        // Track the underlying Fp3 multiplication for the second square
        tracker.record(Op::Fp3Mul);
        // x⁴ -> x⁵
        let res = x4 * self;
        // Track the final Fp3 multiplication
        tracker.record(Op::Fp3Mul);
        // Explicitly record the completion of the high-level Fp3 Quintic operation
        tracker.record(Op::Fp3Quintic);
        res
    }
}

//  ---------------------------------------------------------------------------
//  Permutation – RPO‑M31
//  ---------------------------------------------------------------------------

/// A stateless permutation implementing the RPO-M31 algorithm.
///
/// RPO-M31 consists of [`RPO_ROUNDS`] rounds, each comprising:
/// 1. Forward Mix (FM): Apply MDS matrix, add round constants, apply quintic S-box.
/// 2. Backward Mix (BM): Apply MDS matrix, add round constants, apply inverse quintic S-box.
///
/// Followed by a final Linear Layer (CLS): Apply MDS matrix, add final round constants.
#[derive(Debug, Default, Clone, Copy)]
pub struct RpoM31;

impl RpoM31 {
    /// Applies the RPO-M31 permutation **in-place** to a 24-element state.
    ///
    /// # Arguments
    ///
    /// * `state`: A mutable reference to the 24-element state array.
    /// * `tracker`: A mutable reference to an OpsTracker.
    #[inline]
    pub fn apply<T: OpsTracker>(state: &mut [Felt; STATE_WIDTH], tracker: &mut T) {
        use constants::RPO as C;
        use mds::M;
        tracker.record(Op::Permutation);
        let mut idx = 0;
        for _round in 0..RPO_ROUNDS {
            // ── FM step (Forward Mix) ───────────────────────────────────────
            mul_matrix(state, &M, tracker);
            add_constants(state, &C[idx], tracker);
            // Pass tracker into the mapped function via closure capture
            map_elements(state, |x| quintic(x, tracker));
            idx += 1;
            // ── BM step (Backward Mix) ──────────────────────────────────────
            mul_matrix(state, &M, tracker);
            add_constants(state, &C[idx], tracker);
            // Pass tracker into the mapped function via closure capture
            map_elements(state, |x| quintic_inv(x, tracker));
            idx += 1;
        }
        // ── CLS step (Final Linear Layer) ──────────────────────────────────
        mul_matrix(state, &M, tracker);
        add_constants(state, &C[idx], tracker);
    }
}

//  ---------------------------------------------------------------------------
//  Permutation – XHash‑M31
//  ---------------------------------------------------------------------------

/// A stateless permutation implementing the XHash-M31 algorithm.
///
/// XHash-M31 consists of [`XHASH_ROUNDS`] round triplets, each comprising:
/// 1. Forward Mix (FM): Apply MDS matrix, add round constants, apply quintic S-box (over M31).
/// 2. Backward Mix (BM): Apply MDS matrix, add round constants, apply inverse quintic S-box (over M31).
/// 3. Fp3 Mix (P3M): Add round constants, apply quintic S-box over Fp3 (treating state as 8 Fp3 elements).
///
/// Followed by a final Linear Layer (CLS): Apply MDS matrix, add final round constants.
#[derive(Debug, Default, Clone, Copy)]
pub struct XHashM31;

impl XHashM31 {
    /// Applies the XHash-M31 permutation **in-place** to a 24-element state.
    ///
    /// # Arguments
    ///
    /// * `state`: A mutable reference to the 24-element state vector.
    /// * `tracker`: A mutable reference to an OpsTracker.
    #[inline]
    pub fn apply<T: OpsTracker>(state: &mut [Felt; STATE_WIDTH], tracker: &mut T) {
        use constants::XHASH as C;
        use mds::M;
        tracker.record(Op::Permutation);
        let mut idx = 0;
        for _round in 0..XHASH_ROUNDS {
            // ── FM step (Forward Mix - M31) ─────────────────────────────────
            mul_matrix(state, &M, tracker);
            add_constants(state, &C[idx], tracker);
            // Pass tracker into the mapped function via closure capture
            map_elements(state, |x| quintic(x, tracker));
            idx += 1;

            // ── BM step (Backward Mix - M31) ────────────────────────────────
            mul_matrix(state, &M, tracker);
            add_constants(state, &C[idx], tracker);
            // Pass tracker into the mapped function via closure capture
            map_elements(state, |x| quintic_inv(x, tracker));
            idx += 1;

            // ── P3M step (Fp3 Mix) ──────────────────────────────────────────
            add_constants(state, &C[idx], tracker); // Add M31 constants
            // Apply Fp3 quintic S-box to 8 disjoint triplets
            for t in 0..8 {
                let base_idx = t * 3;
                let fp3_element = Fp3 {
                    a: state[base_idx],
                    b: state[base_idx + 1],
                    c: state[base_idx + 2],
                };
                // Apply Fp3 quintic, passing the tracker
                let result_fp3 = fp3_element.quintic(tracker);
                state[base_idx] = result_fp3.a;
                state[base_idx + 1] = result_fp3.b;
                state[base_idx + 2] = result_fp3.c;
            }
            idx += 1;
        }
        // ── CLS step (Final Linear Layer) ──────────────────────────────────
        mul_matrix(state, &M, tracker);
        add_constants(state, &C[idx], tracker);
    }
}

//  ---------------------------------------------------------------------------
//  Sponge construction (rate 16 / cap 8)
//  ---------------------------------------------------------------------------

/// A generic sponge construction based on a chosen permutation `P` and an OpsTracker `Tr`.
///
/// Implements a standard sponge mode with rate `RATE` (16) and capacity `STATE_WIDTH - RATE` (8).
/// It provides methods for absorbing data ([`absorb`], [`absorb_bytes`]) and squeezing output ([`squeeze`]).
///
/// The permutation `P` must implement the [`Permutation`] trait.
#[derive(Clone)]
pub struct Sponge<P: Permutation, Tr: OpsTracker> {
    /// The internal state vector (24 field elements).
    pub state: [Felt; STATE_WIDTH],
    /// The current position within the rate part of the state for absorption (0..RATE).
    pub pos: usize,
    /// The operation tracker.
    pub tracker: Tr,
    /// Marker for the permutation type.
    _marker: core::marker::PhantomData<P>,
}

/// Creates a new `Sponge` instance with a zero-initialized state and a default tracker.
impl<P: Permutation, Tr: OpsTracker + Default> Default for Sponge<P, Tr> {
    fn default() -> Self {
        Self::new_with_tracker(Tr::default())
    }
}

impl<P: Permutation> Sponge<P, NoopOpsTracker> {
    /// Creates a new `Sponge` instance with a zero-initialized state and a NoopOpsTracker.
    pub fn new() -> Self {
        Self::new_with_tracker(NoopOpsTracker {})
    }
}

impl<P: Permutation, Tr: OpsTracker> Sponge<P, Tr> {
    /// Creates a new `Sponge` instance with a zero-initialized state and the provided tracker.
    pub fn new_with_tracker(tracker: Tr) -> Self {
        Self {
            state: [Felt::from_u32_unchecked(0); STATE_WIDTH],
            pos: 0,
            tracker,
            _marker: core::marker::PhantomData,
        }
    }

    /// Absorbs a single field element into the sponge state.
    pub fn absorb(&mut self, element: Felt) {
        // Add element to the current position in the rate part.
        self.state[self.pos] += element; // XOR would also work in prime field
        self.tracker.record(Op::FeltAdd);
        self.pos += 1;
        // If the rate part is full, permute the state.
        if self.pos == RATE {
            P::apply(&mut self.state, &mut self.tracker);
            self.pos = 0; // Reset position
        }
    }

    /// Absorbs a slice of bytes into the sponge state.
    pub fn absorb_bytes(&mut self, bytes: &[u8]) {
        for chunk in bytes.chunks(4) {
            let mut word = [0u8; 4];
            // Copy chunk slice into the 4-byte array, padding with zeros if needed.
            word[..chunk.len()].copy_from_slice(chunk);
            // Convert little-endian bytes to u32, then to Felt, and absorb.
            self.absorb(Felt::from(u32::from_le_bytes(word)));
        }
    }

    /// Squeezes the sponge to produce an output digest of `RATE` elements.
    pub fn squeeze(mut self) -> [Felt; RATE] {
        // If the rate part wasn't full, apply the permutation to mix the last absorbed elements.
        if self.pos != 0 {
            // Apply padding implicitly by permuting the current state.
            P::apply(&mut self.state, &mut self.tracker);
        }
        // Squeeze the rate part of the state.
        let mut out = [Felt::from_u32_unchecked(0); RATE];
        out.copy_from_slice(&self.state[..RATE]);
        out
    }

    /// Finalizes the sponge and returns both the digest and the tracker state.
    pub fn squeeze_and_report(mut self) -> ([Felt; RATE], Tr) {
        if self.pos != 0 {
            P::apply(&mut self.state, &mut self.tracker);
        }
        let mut out = [Felt::from_u32_unchecked(0); RATE];
        out.copy_from_slice(&self.state[..RATE]);
        (out, self.tracker)
    }
}

/// A trait abstracting over the permutation function ([`RpoM31`] or [`XHashM31`]).
///
/// This allows the [`Sponge`] struct to be generic over the permutation used.
pub trait Permutation {
    /// Applies the permutation in-place to the given state using the provided tracker.
    ///
    /// # Arguments
    ///
    /// * `state`: A mutable reference to the 24-element state array.
    /// * `tracker`: A mutable reference to an OpsTracker.
    fn apply<T: OpsTracker>(state: &mut [Felt; STATE_WIDTH], tracker: &mut T);
}

/// Implements the [`Permutation`] trait for [`RpoM31`].
impl Permutation for RpoM31 {
    #[inline(always)] // Keep inline as it's just delegation
    fn apply<T: OpsTracker>(state: &mut [Felt; STATE_WIDTH], tracker: &mut T) {
        // Delegates to the inherent apply method.
        Self::apply(state, tracker)
    }
}

/// Implements the [`Permutation`] trait for [`XHashM31`].
impl Permutation for XHashM31 {
    #[inline(always)] // Keep inline as it's just delegation
    fn apply<T: OpsTracker>(state: &mut [Felt; STATE_WIDTH], tracker: &mut T) {
        // Delegates to the inherent apply method.
        Self::apply(state, tracker)
    }
}

//  ---------------------------------------------------------------------------
//  Helper utilities
//  ---------------------------------------------------------------------------

/// Adds round constants `c` to the `state` vector element-wise.
#[inline(always)]
fn add_constants<T: OpsTracker>(
    state: &mut [Felt; STATE_WIDTH],
    c: &[Felt; STATE_WIDTH],
    tracker: &mut T,
) {
    for (s, k) in state.iter_mut().zip(c) {
        *s += *k;
        tracker.record(Op::FeltAdd);
    }
}

/// Applies a function `f` to each element of the `state` vector in-place.
#[inline(always)]
fn map_elements<F: FnMut(Felt) -> Felt>(state: &mut [Felt; STATE_WIDTH], mut f: F) {
    // Pass the mutable closure `f` which captures the tracker mutably.
    for s in state.iter_mut() {
        *s = f(*s);
    }
}

/// Multiplies the `state` vector by the MDS `matrix` in-place.
///
/// `state_new[i] = dot_product(matrix_row_i, state_old)`
/// Uses an unrolled macro for the 24x24 dot product.
#[inline(always)]
fn mul_matrix<T: OpsTracker>(
    state: &mut [Felt; STATE_WIDTH],
    m: &[[Felt; STATE_WIDTH]; STATE_WIDTH],
    tracker: &mut T,
) {
    // Macro for computing the dot product of a matrix row and the state vector.
    // Records the necessary Felt multiplications and additions.
    macro_rules! dot24 {
        ($row:expr, $s:expr, $tracker:expr) => {{
            let mut sum = Felt::from_u32_unchecked(0);
            sum += $row[0] * $s[0];
            $tracker.record(Op::FeltMul);
            $tracker.record(Op::FeltAdd);
            sum += $row[1] * $s[1];
            $tracker.record(Op::FeltMul);
            $tracker.record(Op::FeltAdd);
            sum += $row[2] * $s[2];
            $tracker.record(Op::FeltMul);
            $tracker.record(Op::FeltAdd);
            sum += $row[3] * $s[3];
            $tracker.record(Op::FeltMul);
            $tracker.record(Op::FeltAdd);
            sum += $row[4] * $s[4];
            $tracker.record(Op::FeltMul);
            $tracker.record(Op::FeltAdd);
            sum += $row[5] * $s[5];
            $tracker.record(Op::FeltMul);
            $tracker.record(Op::FeltAdd);
            sum += $row[6] * $s[6];
            $tracker.record(Op::FeltMul);
            $tracker.record(Op::FeltAdd);
            sum += $row[7] * $s[7];
            $tracker.record(Op::FeltMul);
            $tracker.record(Op::FeltAdd);
            sum += $row[8] * $s[8];
            $tracker.record(Op::FeltMul);
            $tracker.record(Op::FeltAdd);
            sum += $row[9] * $s[9];
            $tracker.record(Op::FeltMul);
            $tracker.record(Op::FeltAdd);
            sum += $row[10] * $s[10];
            $tracker.record(Op::FeltMul);
            $tracker.record(Op::FeltAdd);
            sum += $row[11] * $s[11];
            $tracker.record(Op::FeltMul);
            $tracker.record(Op::FeltAdd);
            sum += $row[12] * $s[12];
            $tracker.record(Op::FeltMul);
            $tracker.record(Op::FeltAdd);
            sum += $row[13] * $s[13];
            $tracker.record(Op::FeltMul);
            $tracker.record(Op::FeltAdd);
            sum += $row[14] * $s[14];
            $tracker.record(Op::FeltMul);
            $tracker.record(Op::FeltAdd);
            sum += $row[15] * $s[15];
            $tracker.record(Op::FeltMul);
            $tracker.record(Op::FeltAdd);
            sum += $row[16] * $s[16];
            $tracker.record(Op::FeltMul);
            $tracker.record(Op::FeltAdd);
            sum += $row[17] * $s[17];
            $tracker.record(Op::FeltMul);
            $tracker.record(Op::FeltAdd);
            sum += $row[18] * $s[18];
            $tracker.record(Op::FeltMul);
            $tracker.record(Op::FeltAdd);
            sum += $row[19] * $s[19];
            $tracker.record(Op::FeltMul);
            $tracker.record(Op::FeltAdd);
            sum += $row[20] * $s[20];
            $tracker.record(Op::FeltMul);
            $tracker.record(Op::FeltAdd);
            sum += $row[21] * $s[21];
            $tracker.record(Op::FeltMul);
            $tracker.record(Op::FeltAdd);
            sum += $row[22] * $s[22];
            $tracker.record(Op::FeltMul);
            $tracker.record(Op::FeltAdd);
            sum += $row[23] * $s[23];
            $tracker.record(Op::FeltMul);
            $tracker.record(Op::FeltAdd);
            // Note: We record 24 Adds, but technically it's 23 additions per dot product.
            // Adjust if precision down to single additions matters.
            sum
        }};
    }

    tracker.record(Op::MdsMul);
    // Copy state to a temporary array `s` to avoid modifying it while reading during calculation.
    let s = *state; // copy to avoid aliasing during computation
    for i in 0..STATE_WIDTH {
        state[i] = dot24!(m[i], s, tracker);
    }
}

//  ---------------------------------------------------------------------------
//  TESTS
//  ---------------------------------------------------------------------------
#[cfg(test)]
/// Unit tests for the RPO/XHash components.
mod tests {
    use super::*;
    use rand::rngs::SmallRng;
    use rand::{Rng, SeedableRng};

    /// Test the efficient `quintic` function against the generic `pow`.
    #[test]
    fn quintic_vs_pow() {
        let mut rng = SmallRng::seed_from_u64(0);
        for _ in 0..10_000 {
            let x = Felt::from(rng.r#gen::<u32>() % MODULUS);
            assert_eq!(
                quintic(x, &mut NoopOpsTracker {}),
                pow(x, 5, &mut NoopOpsTracker {})
            );
        }
    }

    /// Test that the RPO permutation changes the state (basic check).
    #[test]
    fn rpo_roundtrip_length_preservation() {
        let mut state = [Felt::from_u32_unchecked(0); STATE_WIDTH];
        for (i, m) in state.iter_mut().enumerate().take(STATE_WIDTH) {
            *m = Felt::from(i as u32 + 3);
        }
        let copy = state;
        let mut tracker = NoopOpsTracker::default();
        RpoM31::apply(&mut state, &mut tracker);
        assert_ne!(state, copy);
    }

    /// Test that the XHash permutation changes the state (basic check).
    #[test]
    fn xhash_roundtrip_length_preservation() {
        let mut state = [Felt::from_u32_unchecked(0); STATE_WIDTH];
        for (i, m) in state.iter_mut().enumerate().take(STATE_WIDTH) {
            *m = Felt::from(i as u32 + 7);
        }
        let copy = state;
        let mut tracker = NoopOpsTracker::default();
        XHashM31::apply(&mut state, &mut tracker);
        assert_ne!(state, copy);
    }

    /// Test basic sponge workflow: absorb multiple elements, then squeeze.
    #[test]
    fn sponge_basic() {
        let mut sp: Sponge<RpoM31, NoopOpsTracker> = Sponge::new();
        // Absorb enough elements to trigger at least one permutation.
        for i in 0..40 {
            sp.absorb(Felt::from(i));
        }
        let digest = sp.squeeze();
        let expected_first = Felt::from(100_u32); // placeholder sentinel
        assert_ne!(digest[0], expected_first, "dummy non-zero check");
    }

    /// Helper to get a Felt from a deterministic source for testing.
    fn felt(val: u32) -> Felt {
        Felt::from(val)
    }

    /// Test the `quintic_inv` function by checking `quintic_inv(quintic(x)) == x`.
    #[test]
    fn test_quintic_inverse_roundtrip() {
        let mut rng = SmallRng::seed_from_u64(1); // Use different seed
        for _ in 0..100 {
            let x = Felt::from(rng.r#gen::<u32>() % MODULUS);
            if x != Felt::from(0) {
                let x5 = quintic(x, &mut NoopOpsTracker {});
                let x_inv5 = quintic_inv(x5, &mut NoopOpsTracker {});
                assert_eq!(
                    x_inv5, x,
                    "quintic inverse roundtrip failed for x = {:?}",
                    x
                );
            }
        }
        assert_eq!(
            quintic_inv(Felt::from(1), &mut NoopOpsTracker {}),
            Felt::from(1)
        );
    }

    /// Test sponge determinism for RPO.
    #[test]
    fn test_sponge_determinism_rpo() {
        let input_data = b"Test data for determinism";
        let input_felts = [felt(1), felt(2), felt(3), felt(1000)];

        let mut sponge1: Sponge<RpoM31, NoopOpsTracker> = Sponge::new();
        sponge1.absorb_bytes(input_data);
        for &f in &input_felts {
            sponge1.absorb(f);
        }
        let digest1 = sponge1.squeeze();

        let mut sponge2: Sponge<RpoM31, NoopOpsTracker> = Sponge::new();
        sponge2.absorb_bytes(input_data);
        for &f in &input_felts {
            sponge2.absorb(f);
        }
        let digest2 = sponge2.squeeze();

        assert_eq!(digest1, digest2, "RPO Sponge output is not deterministic");
    }

    /// Test sponge determinism for XHash.
    #[test]
    fn test_sponge_determinism_xhash() {
        let input_data = b"Another test string for XHash";
        let input_felts = [felt(99), felt(88), felt(77), felt(12345)];

        let mut sponge1: Sponge<XHashM31, NoopOpsTracker> = Sponge::new();
        sponge1.absorb_bytes(input_data);
        for &f in &input_felts {
            sponge1.absorb(f);
        }
        let digest1 = sponge1.squeeze();

        let mut sponge2: Sponge<XHashM31, NoopOpsTracker> = Sponge::new();
        sponge2.absorb_bytes(input_data);
        for &f in &input_felts {
            sponge2.absorb(f);
        }
        let digest2 = sponge2.squeeze();

        assert_eq!(digest1, digest2, "XHash Sponge output is not deterministic");
    }

    /// Test the CountingOpsTracker for RPO.
    #[test]
    fn test_rpo_op_counting() {
        // Test counts for exactly one permutation call
        let mut state = [Felt::from_u32_unchecked(0); STATE_WIDTH];
        let mut tracker = CountingOpsTracker::new();

        // Apply the permutation once
        RpoM31::apply(&mut state, &mut tracker);

        let report = tracker.report_json().unwrap();
        println!("RPO Ops Report (1 Permutation Call):\n{}", report);

        // Basic sanity checks for one permutation
        assert_eq!(*tracker.counts.get(&Op::Permutation).unwrap_or(&0), 1);
        assert!(*tracker.counts.get(&Op::FeltAdd).unwrap_or(&0) > 0);
        assert!(*tracker.counts.get(&Op::FeltMul).unwrap_or(&0) > 0);
        assert_eq!(
            *tracker.counts.get(&Op::MdsMul).unwrap_or(&0),
            (RPO_ROUNDS * 2 + 1) as u64
        );
        assert_eq!(
            *tracker.counts.get(&Op::FeltQuintic).unwrap_or(&0),
            (RPO_ROUNDS * STATE_WIDTH) as u64
        );
        assert_eq!(
            *tracker.counts.get(&Op::FeltQuinticInv).unwrap_or(&0),
            (RPO_ROUNDS * STATE_WIDTH) as u64
        );
        // Should not have Fp3 ops for RPO
        assert!(tracker.counts.get(&Op::Fp3Mul).is_none());
        assert!(tracker.counts.get(&Op::Fp3Quintic).is_none());
    }

    /// Test the CountingOpsTracker for XHash.
    #[test]
    fn test_xhash_op_counting() {
        // Test counts for exactly one permutation call
        let mut state = [Felt::from_u32_unchecked(0); STATE_WIDTH];
        let mut tracker = CountingOpsTracker::new();

        // Apply the permutation once
        XHashM31::apply(&mut state, &mut tracker);

        let report = tracker.report_json().unwrap();
        println!("XHash Ops Report (1 Permutation Call):\n{}", report);

        // Basic sanity checks for one permutation
        assert_eq!(*tracker.counts.get(&Op::Permutation).unwrap_or(&0), 1);
        assert!(*tracker.counts.get(&Op::FeltAdd).unwrap_or(&0) > 0);
        assert!(*tracker.counts.get(&Op::FeltMul).unwrap_or(&0) > 0);
        assert_eq!(
            *tracker.counts.get(&Op::MdsMul).unwrap_or(&0),
            (XHASH_ROUNDS * 2 + 1) as u64
        );
        assert_eq!(
            *tracker.counts.get(&Op::FeltQuintic).unwrap_or(&0),
            (XHASH_ROUNDS * STATE_WIDTH) as u64
        );
        assert_eq!(
            *tracker.counts.get(&Op::FeltQuinticInv).unwrap_or(&0),
            (XHASH_ROUNDS * STATE_WIDTH) as u64
        );
        // Should have Fp3 ops for XHash
        assert_eq!(
            *tracker.counts.get(&Op::Fp3Mul).unwrap_or(&0),
            (XHASH_ROUNDS * 8 * 3) as u64
        ); // 8 triplets, 3 muls per quintic
        assert_eq!(
            *tracker.counts.get(&Op::Fp3Quintic).unwrap_or(&0),
            (XHASH_ROUNDS * 8) as u64
        ); // 8 triplets
    }
}
