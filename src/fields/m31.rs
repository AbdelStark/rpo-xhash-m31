use std::fmt::Display;
use std::ops::{
    Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Rem, RemAssign, Sub, SubAssign,
};

use bytemuck::{Pod, Zeroable};
use rand::distributions::{Distribution, Standard};
use serde::{Deserialize, Serialize};

use super::{ComplexConjugate, FieldExpOps};
use crate::impl_field;
pub const MODULUS_BITS: u32 = 31;
pub const N_BYTES_FELT: usize = 4;
pub const P: u32 = 2147483647; // 2 ** 31 - 1

#[repr(transparent)]
#[derive(
    Copy, Clone, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize,
)]
pub struct M31(pub u32);
pub type BaseField = M31;

unsafe impl Pod for M31 {}

unsafe impl Zeroable for M31 {
    fn zeroed() -> Self {
        unsafe { core::mem::zeroed() }
    }
}

impl_field!(M31, P);

impl M31 {
    /// Returns `val % P` when `val` is in the range `[0, 2P)`.
    pub fn partial_reduce(val: u32) -> Self {
        Self(val.checked_sub(P).unwrap_or(val))
    }

    /// Returns `val % P` when `val` is in the range `[0, P^2)`.
    pub const fn reduce(val: u64) -> Self {
        Self((((((val >> MODULUS_BITS) + val + 1) >> MODULUS_BITS) + val) & (P as u64)) as u32)
    }

    pub const fn from_u32_unchecked(arg: u32) -> Self {
        Self(arg)
    }

    pub fn inverse(&self) -> Self {
        assert!(!self.is_zero(), "0 has no inverse");
        pow2147483645(*self)
    }
}

impl Display for M31 {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.0)
    }
}

impl Add for M31 {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self::partial_reduce(self.0 + rhs.0)
    }
}

impl Neg for M31 {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self::partial_reduce(P - self.0)
    }
}

impl Sub for M31 {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        Self::partial_reduce(self.0 + P - rhs.0)
    }
}

impl Mul for M31 {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        Self::reduce((self.0 as u64) * (rhs.0 as u64))
    }
}

impl FieldExpOps for M31 {
    fn inverse(&self) -> Self {
        self.inverse()
    }
}

impl ComplexConjugate for M31 {
    fn complex_conjugate(&self) -> Self {
        *self
    }
}

impl One for M31 {
    fn one() -> Self {
        Self(1)
    }
}

impl Zero for M31 {
    fn zero() -> Self {
        Self(0)
    }

    fn is_zero(&self) -> bool {
        *self == Self::zero()
    }
}

impl From<usize> for M31 {
    fn from(value: usize) -> Self {
        M31::reduce(value.try_into().unwrap())
    }
}

impl From<u32> for M31 {
    fn from(value: u32) -> Self {
        M31::reduce(value.into())
    }
}

impl From<i32> for M31 {
    fn from(value: i32) -> Self {
        if value < 0 {
            const P2: u64 = 2 * P as u64;
            return M31::reduce(P2 - value.unsigned_abs() as u64);
        }

        M31::reduce(value.unsigned_abs() as u64)
    }
}

impl Distribution<M31> for Standard {
    // Not intended for cryptographic use. Should only be used in tests and benchmarks.
    fn sample<R: rand::Rng + ?Sized>(&self, rng: &mut R) -> M31 {
        M31(rng.gen_range(0..P))
    }
}

#[cfg(test)]
#[macro_export]
macro_rules! m31 {
    ($m:expr) => {
        $crate::fields::m31::M31::from_u32_unchecked($m)
    };
}

/// Computes `v^((2^31-1)-2)`.
///
/// Computes the multiplicative inverse of [`M31`] elements with 37 multiplications vs naive 60
/// multiplications. Made generic to support both vectorized and non-vectorized implementations.
/// Multiplication tree found with [addchain](https://github.com/mmcloughlin/addchain).
///
pub fn pow2147483645<T: FieldExpOps>(v: T) -> T {
    let t0 = sqn::<2, T>(v.clone()) * v.clone();
    let t1 = sqn::<1, T>(t0.clone()) * t0.clone();
    let t2 = sqn::<3, T>(t1.clone()) * t0.clone();
    let t3 = sqn::<1, T>(t2.clone()) * t0.clone();
    let t4 = sqn::<8, T>(t3.clone()) * t3.clone();
    let t5 = sqn::<8, T>(t4.clone()) * t3.clone();
    sqn::<7, T>(t5) * t2
}

/// Computes `v^(2*n)`.
fn sqn<const N: usize, T: FieldExpOps>(mut v: T) -> T {
    for _ in 0..N {
        v = v.square();
    }
    v
}
