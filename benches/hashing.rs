use criterion::{Criterion, black_box, criterion_group, criterion_main};
use rand::{RngCore, SeedableRng, rngs::SmallRng};
use rpo_xhash_m31::{Felt, RpoM31, STATE_WIDTH, Sponge, XHashM31};

fn create_random_state(rng: &mut SmallRng) -> [Felt; STATE_WIDTH] {
    let mut state = [Felt::default(); STATE_WIDTH];
    for elem in state.iter_mut() {
        *elem = Felt::from(rng.next_u32());
    }
    state
}

fn bench_permutations(c: &mut Criterion) {
    let mut group = c.benchmark_group("Permutations");
    let mut rng = SmallRng::seed_from_u64(42);

    group.bench_function("RpoM31::apply", |b| {
        let mut state = create_random_state(&mut rng);
        b.iter(|| RpoM31::apply(black_box(&mut state)))
    });

    group.bench_function("XHashM31::apply", |b| {
        let mut state = create_random_state(&mut rng);
        b.iter(|| XHashM31::apply(black_box(&mut state)))
    });

    group.finish();
}

fn bench_sponge(c: &mut Criterion) {
    let mut group = c.benchmark_group("Sponge Operations");
    let mut rng = SmallRng::seed_from_u64(42);
    let data_sizes = [64, 256, 1024, 4096]; // Bytes

    for size in data_sizes.iter() {
        let mut input_data = vec![0u8; *size];
        rng.fill_bytes(&mut input_data);

        group.bench_with_input(
            criterion::BenchmarkId::new("RpoM31::absorb_bytes", size),
            &input_data,
            |b, data| {
                b.iter(|| {
                    let mut sponge = Sponge::<RpoM31>::new();
                    sponge.absorb_bytes(black_box(data));
                    black_box(sponge); // Prevent elimination
                })
            },
        );

        group.bench_with_input(
            criterion::BenchmarkId::new("XHashM31::absorb_bytes", size),
            &input_data,
            |b, data| {
                b.iter(|| {
                    let mut sponge = Sponge::<XHashM31>::new();
                    sponge.absorb_bytes(black_box(data));
                    black_box(sponge); // Prevent elimination
                })
            },
        );
    }

    let mut rpo_sponge_filled = Sponge::<RpoM31>::new();
    let mut xh_sponge_filled = Sponge::<XHashM31>::new();
    let mut large_input = vec![0u8; 1024]; // Arbitrary data to fill sponge
    rng.fill_bytes(&mut large_input);
    rpo_sponge_filled.absorb_bytes(&large_input);
    xh_sponge_filled.absorb_bytes(&large_input);

    group.bench_function("RpoM31::squeeze", |b| {
        // Clone sponge inside iter to reset state for each measurement
        b.iter_batched(
            || rpo_sponge_filled.clone(),
            |sponge| black_box(sponge.squeeze()),
            criterion::BatchSize::SmallInput,
        )
    });

    group.bench_function("XHashM31::squeeze", |b| {
        b.iter_batched(
            || xh_sponge_filled.clone(),
            |sponge| black_box(sponge.squeeze()),
            criterion::BatchSize::SmallInput,
        )
    });

    group.finish();
}

criterion_group!(benches, bench_permutations, bench_sponge);
criterion_main!(benches);
