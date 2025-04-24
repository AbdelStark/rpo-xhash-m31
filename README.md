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
Zero Knowledge proof systems.  *Arithmetisation‑Oriented* (AO) primitives fix this by favouring low‑degree, highly
regular operations.

* **RPO‑M31** – a Rescue‑Prime Optimised permutation adapted to the 31‑bit Mersenne field.
* **XHash‑M31** – interleaves RPO rounds with cubic‑extension S‑box layers for extra diffusion.

Both operate on a 24‑element state → 16‑element rate / 8‑element capacity, yielding **~124‑bit**
generic security (Section 3 of the paper).

---

## 🎮 Quick‑start

```bash
# add the dependency (until published use a git URL)
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
# 1. run unit + integration tests
cargo test --all-targets --release

# 2. benchmark (criterion)
cargo bench
```

---

## 📄 License

Licensed under the [MIT license](LICENSE).  

---

## 📖 References

* [RPO-M31 and XHash-M31: Efficient Hash Functions for Circle STARKs](https://eprint.iacr.org/2024/1635.pdf)
* [Circle STARKs paper](https://eprint.iacr.org/2024/278)
* [STWO prover](https://github.com/starkware-libs/stwo)

---

> “Simplicity is prerequisite for reliability.” ― **E. W. Dijkstra**
