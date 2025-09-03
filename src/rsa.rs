use crate::util::*;
use num_bigint::{BigInt, BigUint, RandBigInt, ToBigInt};
use num_traits::{One, Zero};

// OpenSSL style names for the struct fields
pub struct PublicKey {
    pub n: BigUint, // modulus
    pub e: BigUint, // public exponent
}

pub struct PrivateKey {
    pub n: BigUint, // modulus
    pub e: BigUint, // public exponent
    pub d: BigUint, // private exponent
    pub p: BigUint, // prime factor
    pub q: BigUint, // prime factor
    pub dmp: BigUint, // d mod (p-1)
    pub dmq: BigUint, // d mod (q-1)
    pub iqmp: BigUint, // q^(-1) mod p
}

#[derive(Debug)]
pub enum RsaError {
    MessageOutOfRange,
    CiphertextOutOfRange,
    NonInvertibleExponent,
    PrimeGenFailed,
    InvalidBitLength,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::util::*;
    use num_bigint::{BigInt, BigUint, ToBigInt, ToBigUint};
    use num_traits::{One, Zero};
    use num_bigint::RandBigInt;
    use rand_chacha::{ChaCha20Rng, rand_core::SeedableRng};

    // helper functions
}
