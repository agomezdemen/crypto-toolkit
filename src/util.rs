// Alejandro Gomez De Mendieta util.rs
// Utility functions for number theory and cryptographic math
// Includes modular exponentiation, GCD, and modular inverse

// num-bigint crate for RSA-scale numbers
use num_bigint::{BigInt, BigUint, RandBigInt, ToBigInt};
use num_integer::Integer;
use num_traits::{One, Zero};


// Performs modular exponentiation
pub fn modular_expo(base: &BigUint, exponent: &BigUint, modulus: &BigUint) -> BigUint{
    let result: BigUint = 0;
    
    // covering modulus = 1 edge case for efficiency
    if modulus == 1 {
        return result;
    }
    


}
