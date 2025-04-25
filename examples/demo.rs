//! Simple CLI demo.
//!
//! `$ cargo run --release --example demo "starks"`
//!
//! Prints two 64-byte (128 field-bit) digests – one produced with **RPO-M31**,
//! the other with **XHash-M31**.

use rpo_xhash_m31::{Felt, NoopOpsTracker, RpoM31, Sponge, XHashM31};
use std::env;

/// Helper – convert the 16-element digest into a lowercase hex string.
fn digest_to_hex(words: &[Felt; 16]) -> String {
    let mut out = String::with_capacity(16 * 8);
    for w in words {
        out.push_str(&format!("{:08x}", w.0)); // M31 fits in 32 bits.
    }
    out
}

fn main() {
    // ----------------------------------------------- read CLI argument
    let input = env::args().nth(1).expect("Please provide a message");
    let bytes = input.as_bytes();

    // ----------------------------------------------- RPO-M31
    let mut rpo = Sponge::<RpoM31, NoopOpsTracker>::new();
    rpo.absorb_bytes(bytes);
    let rpo_digest = rpo.squeeze();
    println!("RPO-M31  : {}", digest_to_hex(&rpo_digest));

    // ----------------------------------------------- XHash-M31
    let mut xh = Sponge::<XHashM31, NoopOpsTracker>::new();
    xh.absorb_bytes(bytes);
    let xh_digest = xh.squeeze();
    println!("XHash-M31: {}", digest_to_hex(&xh_digest));
}
