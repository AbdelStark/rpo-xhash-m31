use rpo_xhash_m31::{CountingOpsTracker, Felt, Op, RpoM31, STATE_WIDTH, XHashM31};
use std::collections::HashMap;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // vByte cost estimates per operation (as defined in the plan)
    let costs: HashMap<Op, u64> = [
        (Op::FeltAdd, 10),
        (Op::FeltMul, 400),
        (Op::FeltSquare, 300),       // Optimized Mul
        (Op::FeltQuintic, 1000),     // 2*Square + 1*Mul = 2*300 + 400
        (Op::FeltQuinticInv, 15700), // Approx 31*Square + 16*Mul = 31*300 + 16*400
        (Op::MdsMul, 236160),        // 576*Mul + 576*Add = 576*400 + 576*10
        (Op::Fp3Mul, 3350),          // Approx 8*Mul + 15*Add = 8*400 + 15*10
        (Op::Fp3Quintic, 10050),     // 3*Fp3Mul = 3*3350
    ]
    .iter()
    .cloned()
    .collect();

    // --- Dynamically Calculate RPO Counts ---
    let mut rpo_state = [Felt::from_u32_unchecked(0); STATE_WIDTH];
    let mut rpo_tracker = CountingOpsTracker::new();
    RpoM31::apply(&mut rpo_state, &mut rpo_tracker);
    let rpo_counts = rpo_tracker.get_counts(); // Use the getter method

    // --- Calculate RPO Cost ---
    let mut total_rpo_cost: u64 = 0;
    println!("--- RPO Cost Breakdown ---");
    // Sort for consistent output order
    let mut rpo_ops: Vec<_> = rpo_counts.keys().collect();
    rpo_ops.sort_by_key(|op| format!("{:?}", op));

    for op in rpo_ops {
        let count = rpo_counts.get(op).unwrap(); // We know the key exists
        if let Some(cost) = costs.get(op) {
            let op_total_cost = count * cost;
            total_rpo_cost += op_total_cost;
            println!("{:?}: {} * {} = {} vbytes", op, count, cost, op_total_cost);
        } else if *op != Op::Permutation {
            // Ignore Permutation op itself
            println!("Warning: Cost not defined for RPO Op: {:?}", op);
        }
    }
    println!("--------------------------");
    println!(
        "Total Estimated RPO Permutation Cost: {} vbytes",
        total_rpo_cost
    );
    println!(
        "                                      ~{} KB",
        total_rpo_cost / 1024
    );

    // --- Dynamically Calculate XHash Counts ---
    let mut xhash_state = [Felt::from_u32_unchecked(0); STATE_WIDTH];
    let mut xhash_tracker = CountingOpsTracker::new();
    XHashM31::apply(&mut xhash_state, &mut xhash_tracker);
    let xhash_counts = xhash_tracker.get_counts(); // Use the getter method

    // --- Calculate XHash Cost ---
    let mut total_xhash_cost: u64 = 0;
    println!("\n--- XHash Cost Breakdown ---");
    // Sort for consistent output order
    let mut xhash_ops: Vec<_> = xhash_counts.keys().collect();
    xhash_ops.sort_by_key(|op| format!("{:?}", op));

    for op in xhash_ops {
        let count = xhash_counts.get(op).unwrap(); // We know the key exists
        if let Some(cost) = costs.get(op) {
            let op_total_cost = count * cost;
            total_xhash_cost += op_total_cost;
            println!("{:?}: {} * {} = {} vbytes", op, count, cost, op_total_cost);
        } else if *op != Op::Permutation {
            // Ignore Permutation op itself
            println!("Warning: Cost not defined for XHash Op: {:?}", op);
        }
    }
    println!("--------------------------");
    println!(
        "Total Estimated XHash Permutation Cost: {} vbytes",
        total_xhash_cost
    );
    println!(
        "                                       ~{} KB",
        total_xhash_cost / 1024
    );

    Ok(())
}
