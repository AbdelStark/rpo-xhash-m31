#!/usr/bin/env rust-script
// Dependencies can be added inline:
// cargo-deps: serde="1.0", serde_json="1.0"

use serde::{Deserialize, Serialize};
use std::collections::HashMap;

// Re-define the Op enum for parsing the report JSON
// (Alternatively, could restructure the main lib to share this)
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
enum Op {
    FeltAdd,
    FeltMul,
    FeltSquare,
    FeltQuintic,
    FeltQuinticInv,
    MdsMul,
    Fp3Add, // Keep for potential future refinement
    Fp3Sub, // Keep for potential future refinement
    Fp3Mul,
    Fp3Quintic,
    Permutation, // Ignored in cost calculation
}

// Operation counts obtained from `cargo test -- --nocapture test_..._op_counting`
const RPO_COUNTS_JSON: &str = r#"
{
  "FeltAdd": 9000,
  "MdsMul": 15,
  "FeltSquare": 5544,
  "FeltQuintic": 168,
  "Permutation": 1,
  "FeltQuinticInv": 168,
  "FeltMul": 11496
}
"#;

const XHASH_COUNTS_JSON: &str = r#"
{
  "FeltQuinticInv": 72,
  "FeltMul": 5256,
  "FeltAdd": 4272,
  "Fp3Quintic": 24,
  "MdsMul": 7,
  "Fp3Mul": 72,
  "FeltSquare": 2376,
  "Permutation": 1,
  "FeltQuintic": 72
}
"#;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // vByte cost estimates per operation (as defined in the plan)
    let costs: HashMap<Op, u64> = [
        (Op::FeltAdd, 10),
        (Op::FeltMul, 400),
        (Op::FeltSquare, 300), // Optimized Mul
        // High-level op costs derived from base ops:
        (Op::FeltQuintic, 1000),     // 2*Square + 1*Mul = 2*300 + 400
        (Op::FeltQuinticInv, 15700), // Approx 31*Square + 16*Mul = 31*300 + 16*400
        (Op::MdsMul, 236160),        // 576*Mul + 576*Add = 576*400 + 576*10
        (Op::Fp3Mul, 3350),          // Approx 8*Mul + 15*Add = 8*400 + 15*10
        (Op::Fp3Quintic, 10050),     // 3*Fp3Mul = 3*3350
    ]
    .iter()
    .cloned()
    .collect();

    // --- Calculate RPO Cost ---
    let rpo_counts: HashMap<Op, u64> = serde_json::from_str(RPO_COUNTS_JSON)?;
    let mut total_rpo_cost: u64 = 0;

    println!("--- RPO Cost Breakdown ---");
    for (op, count) in &rpo_counts {
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

    // --- Calculate XHash Cost ---
    let xhash_counts: HashMap<Op, u64> = serde_json::from_str(XHASH_COUNTS_JSON)?;
    let mut total_xhash_cost: u64 = 0;

    println!(
        "
--- XHash Cost Breakdown ---"
    );
    for (op, count) in &xhash_counts {
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
