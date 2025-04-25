# rpo-xhash-m31

> A Rust implementation of the RPO‑M31 and XHash‑M31 hash functions for Circle‑STARKs.
>
> Paper: <https://eprint.iacr.org/2024/1635.pdf>

[![Crates.io](https://img.shields.io/crates/v/rpo-xhash-m31.svg)](https://crates.io/crates/rpo-xhash-m31)
[![Docs.rs](https://docs.rs/rpo-xhash-m31/badge.svg)](https://docs.rs/rpo-xhash-m31)
![License MIT](https://img.shields.io/badge/License-MIT-green)

---

## 📚 Background

Traditional hash functions such as SHA‑2/SHA‑3 are **algebraically complex** and inefficient in
Zero Knowledge proof systems. _Arithmetisation‑Oriented_ (AO) primitives fix this by favouring low‑degree, highly
regular operations.

- **RPO‑M31** – a Rescue‑Prime Optimised permutation adapted to the 31‑bit Mersenne field.
- **XHash‑M31** – interleaves RPO rounds with cubic‑extension S‑box layers for extra diffusion.

Both operate on a 24‑element state → 16‑element rate / 8‑element capacity, yielding **~124‑bit**
generic security (Section 3 of the paper).

---

## 🎮 Quick‑start

```bash
# add the dependency
cargo add rpo-xhash-m31

# run the demo
cargo run --release --example demo "starks"
```

Sample output:

```text
RPO-M31  : 7eddbb721cfcc3ea2401dc0601b966783f350d6814b46f71513485970df4ec80614d3a3c1537f0262c5c839d05511b011c2b196611613b80383cbd127b95b2a3
XHash-M31: 4904db4703a3417e01120c542bee834f1d9c96c005dc14065d55234a234885ac6708a188495b831c4eb2732c73c886392ff6d95660dced5b26d598bd7c13f879
```

### Library usage

```rust
// Import the library elements
use rpo_xhash_m31::{Sponge, RpoM31, XHashM31};

// Create an input
let input = "starks".to_string();
let bytes = input.as_bytes();

// ----------------------------------------------- RPO-M31
let mut rpo = Sponge::<RpoM31>::new();
rpo.absorb_bytes(bytes);
let rpo_digest = rpo.squeeze();

// ----------------------------------------------- XHash-M31
let mut xh = Sponge::<XHashM31>::new();
xh.absorb_bytes(bytes);
let xh_digest = xh.squeeze();
```

## 🛠️ Development

```bash
# 1. run unit + integration tests
cargo test --all-targets --release

# 2. benchmark (criterion)
cargo bench
```

---

## 💡 Bitcoin Script Cost Estimation

To evaluate the feasibility of implementing these hash functions directly in Bitcoin Script, a preliminary cost estimation in virtual bytes (vbytes) was performed.

### Methodology

1. **Operation Counting:** The core permutation logic (`RpoM31::apply` and `XHashM31::apply`) was instrumented using an `OpsTracker` to count the exact number of low-level field operations executed during a single permutation.
2. **vByte Cost Assignment:** Each tracked operation was assigned an estimated vByte cost. These estimates are based on the assumption of using lookup tables for M31 field arithmetic within Bitcoin Script. Key estimates include:
    - `FeltAdd`: 10 vbytes
    - `FeltMul`: 400 vbytes
    - `FeltSquare`: 300 vbytes (Optimized multiplication)
    - `FeltQuintic`: 1000 vbytes (Derived: 2*Square + 1*Mul)
    - `FeltQuinticInv`: 15700 vbytes (Derived: Based on `pow` logic, ~31*Square + ~16*Mul)
    - `MdsMul`: 236160 vbytes (Derived: 576*Mul + 576*Add)
    - `Fp3Mul`: 3350 vbytes (Derived: Approx 8*Mul + 15*Add)
    - `Fp3Quintic`: 10050 vbytes (Derived: 3\*Fp3Mul)
3. **Total Cost Calculation:** The total estimated cost for one permutation was calculated by summing the product of the operation count and its estimated vByte cost (`total_cost = Σ (count[op] * cost[op])`).

The exact counts and the calculation script can be found in `src/bin/cost_estimator.rs`.

### Estimated Costs (Per Permutation)

- **RPO-M31:** ~12,700,000 vbytes (~12.4 MB)
- **XHash-M31:** ~6,200,000 vbytes (~6.1 MB)

### Conclusion

The estimated vByte costs are exceedingly high, primarily due to the extensive use of M31 multiplications within the MDS matrix layer and the S-box inversions. Implementing either RPO-M31 or XHash-M31 directly in Bitcoin Script using current known techniques (like lookup tables) appears infeasible as the cost dramatically exceeds practical transaction and block size limits.

---

## 📄 License

Licensed under the [MIT license](LICENSE).

---

## 📖 References

- [RPO-M31 and XHash-M31: Efficient Hash Functions for Circle STARKs](https://eprint.iacr.org/2024/1635.pdf)
- [Circle STARKs paper](https://eprint.iacr.org/2024/278)
- [STWO prover](https://github.com/starkware-libs/stwo)

---

> "Simplicity is prerequisite for reliability." ― **E. W. Dijkstra**
