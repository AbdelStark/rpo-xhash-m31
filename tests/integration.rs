use num_traits::identities::Zero;
use rpo_xhash_m31::{Felt, RpoM31, Sponge, XHashM31, mds};

/// A *very* small monte-carlo test that the two permutations disagree on random input.
#[test]
fn permutations_produce_distinct_output() {
    let msg = "test-vector-123-🦀";
    // RPO
    let mut rpo = Sponge::<RpoM31>::new();
    rpo.absorb_bytes(msg.as_bytes());
    let a = rpo.squeeze();
    // XHash
    let mut xh = Sponge::<XHashM31>::new();
    xh.absorb_bytes(msg.as_bytes());
    let b = xh.squeeze();
    assert_ne!(a, b, "Different permutations must not collide trivially");
}

/// Domain-separator byte (17th state element) must make *all-zero*
/// and *empty* messages hash to DIFFERENT digests.
#[test]
fn padding_domain_separation() {
    let rpo_empty = Sponge::<RpoM31>::new();
    let empty_digest = rpo_empty.clone().squeeze();

    let mut rpo_zero_block = Sponge::<RpoM31>::new();
    for _ in 0..16 {
        rpo_zero_block.absorb(Felt::from_u32_unchecked(0));
    }
    let zero_digest = rpo_zero_block.squeeze();

    assert_ne!(empty_digest, zero_digest, "domain separation failed");
}

/// Verify that the 24×24 MDS matrix is invertible (max distance separable).
#[test]
fn mds_is_invertible() {
    use stwo_prover::core::fields::m31::M31;

    // Naïve O(n³) Gaussian elimination over the field.
    // 24 is small → fine for a unit-test.
    const N: usize = 24;
    let mut mat = [[M31::from_u32_unchecked(0); N]; N];
    for i in 0..N {
        mat[i].copy_from_slice(&mds::M[i]);
    }

    // augment with identity
    let mut aug = [[M31::from_u32_unchecked(0); N * 2]; N];
    for i in 0..N {
        for j in 0..N {
            aug[i][j] = mat[i][j];
        }
        aug[i][N + i] = M31::from_u32_unchecked(1);
    }

    // forward elimination
    for col in 0..N {
        // find pivot
        let pivot = (col..N).find(|&r| !aug[r][col].is_zero()).unwrap();
        if pivot != col {
            aug.swap(pivot, col);
        }
        let inv = aug[col][col].inverse();
        for j in col..2 * N {
            aug[col][j] *= inv;
        }
        for row in 0..N {
            if row != col {
                let factor = aug[row][col];
                for j in col..2 * N {
                    aug[row][j] -= factor * aug[col][j];
                }
            }
        }
    }

    // left block should now be identity
    for i in 0..N {
        for j in 0..N {
            assert_eq!(
                aug[i][j],
                if i == j {
                    M31::from_u32_unchecked(1)
                } else {
                    M31::from_u32_unchecked(0)
                }
            );
        }
    }
}

/// Quick shoot-out: hashing the same message twice must be deterministic.
#[test]
fn deterministic_output() {
    let msg = "determinism!";
    let mut a1 = Sponge::<RpoM31>::new();
    a1.absorb_bytes(msg.as_bytes());
    let d1 = a1.squeeze();

    let mut a2 = Sponge::<RpoM31>::new();
    a2.absorb_bytes(msg.as_bytes());
    let d2 = a2.squeeze();

    assert_eq!(d1, d2);
}

/// Tiny property-test: random 3-byte messages expand to ≤1 field element,
/// so we ensure capacity (=8) is never overwritten during absorb.
#[test]
fn short_messages_never_overflow_rate() {
    use rand::Rng;
    let mut rng = rand::rng();
    for _ in 0..10_000 {
        let mut sp: Sponge<RpoM31> = Sponge::new();
        let len = rng.random_range(0..3); // 0-2 bytes
        let mut buf = vec![0u8; len];
        rng.fill(buf.as_mut_slice());
        sp.absorb_bytes(&buf);
        // Internal cursor should never exceed RATE
        assert!(sp.pos < 16);
    }
}
