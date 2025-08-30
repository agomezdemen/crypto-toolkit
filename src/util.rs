// Alejandro Gomez De Mendieta util.rs
// Utility functions for number theory and cryptographic math
// Includes modular exponentiation, GCD, and modular inverse

use std::fmt;
// num-bigint crate for RSA-scale numbers
use num_bigint::{BigInt, BigUint, RandBigInt, ToBigInt};
use num_integer::Integer;
use num_traits::{One, Zero};

// special error type for future extensibility and idiomatic rust conventions
#[derive(Debug)]
pub enum ModExpError {
    InvalidModulus, // for zero case
}
// implementing Display for nice error printing
impl fmt::Display for ModExpError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // match for future extensibility if necessary
        match self {
            ModExpError::InvalidModulus => write!(f, "modulus must be greater than 0"),
        }
    }
}
// implementing Error to integrate with ? operator and crates like anyhow
impl std::error::Error for ModExpError {}

// Performs modular exponentiation using right-to-left binary method
// this method of modular exponentiation reduces time complexity from O(exponent) (naive)
// to O(log(exponent)) which is necessary when working with large numbers required by RSA
pub fn modexp(
    base: &BigUint,
    exponent: &BigUint,
    modulus: &BigUint,
) -> Result<BigUint, ModExpError> {
    // modulus = 0 is an invalid case
    if modulus.is_zero() {
        return Err(ModExpError::InvalidModulus);
    }

    // covering modulus = 1 edge case
    if *modulus == BigUint::one() {
        return Ok(BigUint::ZERO);
    }

    // covering exponent = 0 edge case
    if *exponent == BigUint::ZERO {
        return Ok(BigUint::one());
    }

    let mut result = BigUint::one();
    let mut expo_prime = exponent.clone();

    // the base is reduced here to make computations easier because
    // congruent results repeat when the base is bigger than the modulus
    let mut b = base % modulus;

    // exponentiation process begins
    while !expo_prime.is_zero() {
        // if the exponent is odd (LSB = 1): a^(2q+1) = (a^2q) * a
        if expo_prime.is_odd() {
            result = (&result * &b) % modulus;
        }

        // dropping the lowest bit of the exponent (divides by 2)
        // halves the problem size each step which is what gives us O(log(exponent))
        expo_prime >>= 1;

        // the base at this step represents a^(2^t) for the next bit
        // allowing us to simply square the base
        b = (&b * &b) % modulus;
    }

    Ok(result)
}

// Implements binary GCD algorithm (Steins algorithm)
// for maximum performance to find greatest common divisor
pub fn gcd(a: &BigUint, b: &BigUint) -> BigUint {
    let mut u = a.clone();
    let mut v = b.clone();
    
    // checking edge case where if either is zero then the GCD is the non-zero variable
    if u.is_zero() {
        return v;
    }
    if v.is_zero() {
        return u;
    }
    
    let mut k: usize = 0;

    // this loop counts the common power of two in variable u and v
    // this loop also removes common factors of two because that will not affect the GCD
    // and will help reduce the computations needed
    while u.is_even() && v.is_even() {
        u >>= 1;
        v >>= 1;
        k += 1;
    }

    // making variable u into an odd number to
    // remove remaining factors of two from variable u
    while u.is_even() {
        u >>= 1;
    }

    // this loop is where the main reduction happens
    while !v.is_zero() {

        // making variable v into an odd number to
        // remove remaining factors of two from variable v
        while v.is_even() {
            v >>= 1;
        }
        
        // if statement to swap variable so that the subtraction is always positive
        if u > v {
            std::mem::swap(&mut u, &mut v)
        }
        
        // here we take advantage of the fact 
        // that gcd(u, v) = gcd(u, v - u) when v >= u
        v = &v - &u;
    }

    // here the common power which was factored
    // from the first while loop is restored
    // after all the reductions and the restored common
    // factor u is the GCD for the odd reduced pair
    u << k
}

// Performs the extended greatest common divisor process
pub fn egcd() {}

// Performs modular inverse process
pub fn modinv() {}

#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::BigUint;

    #[test]
    fn modexp_basic_examples() {
        // 5^117 % 97 = 36
        let b = BigUint::from(5u32);
        let e = BigUint::from(117u32);
        let m = BigUint::from(97u32);
        let result = modexp(&b, &e, &m).expect("valid modulus");
        assert_eq!(result, BigUint::from(77u32));

        // 2^10 % 1000 = 24
        let b = BigUint::from(2u32);
        let e = BigUint::from(10u32);
        let m = BigUint::from(1000u32);
        let result = modexp(&b, &e, &m).expect("valid modulus");
        assert_eq!(result, BigUint::from(24u32));
    }

    #[test]
    fn modexp_zero_exponent() {
        // b^0 % m = 1 for any b, m>0
        let b = BigUint::from(7u32);
        let e = BigUint::ZERO;
        let m = BigUint::from(13u32);
        let result = modexp(&b, &e, &m).expect("valid modulus");
        assert_eq!(result, BigUint::one());
    }

    #[test]
    fn modexp_modulus_one() {
        // anything % 1 = 0
        let b = BigUint::from(123u32);
        let e = BigUint::from(456u32);
        let m = BigUint::one();
        let result = modexp(&b, &e, &m).expect("valid modulus");
        assert_eq!(result, BigUint::ZERO);
    }

    #[test]
    fn modexp_invalid_modulus() {
        // result should be InvalidModulus error
        // when modulus = 0
        let b = BigUint::from(3u32);
        let e = BigUint::from(5u32);
        let m = BigUint::ZERO;
        let error = modexp(&b, &e, &m).unwrap_err();
        assert_eq!(error.to_string(), "modulus must be greater than 0");
    }

    #[test]
    fn modexp_modpow_match_with_randoms() {
        use num_bigint::RandBigInt;
        use rand_chacha::{ChaCha20Rng, rand_core::SeedableRng};

        let mut rng = ChaCha20Rng::from_seed([42u8; 32]);

        // random trials with 256-bit values
        for _ in 0..5 {
            let b = rng.gen_biguint(256);
            let e = rng.gen_biguint(256);
            let mut m = rng.gen_biguint(256);
            // ensure modulus is greater than 1
            if m <= BigUint::one() {
                m = BigUint::from(2u32);
            }

            let my_result = modexp(&b, &e, &m).expect("valid modulus");
            let reference = b.modpow(&e, &m);
            assert_eq!(my_result, reference);
        }
    }


}
