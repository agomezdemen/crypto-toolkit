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

// Performs the iterative extended euclidean algorithm
// returns (gcd, x, y) which satisfy a*x + b*y = gcd
pub fn egcd(a: &BigInt, b: &BigInt) -> (BigInt, BigInt, BigInt) {
    let mut old_r = a.clone();
    let mut r = b.clone();

    let mut old_s = BigInt::one();
    let mut s = BigInt::ZERO;

    let mut old_t = BigInt::ZERO;
    let mut t = BigInt::one();

    while !r.is_zero() {
        // old_r = q*r + remainder
        let (q, remainder) = old_r.div_rem(&r);
        old_r = r;
        r = remainder;

        let temp_s = old_s - &q * &s;
        old_s = s;
        s = temp_s;

        let temp_t = old_t - &q * &t;
        old_t = t;
        t = temp_t;
    }

    (old_r, old_s, old_t)
}

// Performs modular inverse
pub fn modinv(a: &BigInt, m: &BigInt) -> Option<BigInt> {
    // egcd takes care of the most important calculations
    // needed to determine if there is an inverse
    let (gcd, mut x, _y) = egcd(a, m);

    if gcd == -(BigInt::one()) {
        x = -(x);
    }
    // if the gcd is not one then there is no inverse
    else if gcd != BigInt::one() {
        return None;
    }

    // normalizing keeps results betewen [0, m)
    // otherwise we would get negative ansers
    let normalized_result = ((x % m) + m) % m;
    Some(normalized_result)
}

// ---------------------tests--------------------------
#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::BigUint;

    mod modexp_tests {
        use super::*;

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

    mod gcd_tests {
        use super::*;

        #[test]
        fn gcd_zero_cases() {
            assert_eq!(gcd(&BigUint::ZERO, &BigUint::ZERO), BigUint::ZERO);
            assert_eq!(
                gcd(&BigUint::ZERO, &BigUint::from(58u32)),
                BigUint::from(58u32)
            );
            assert_eq!(
                gcd(&BigUint::from(23u32), &BigUint::ZERO),
                BigUint::from(23u32)
            );
        }

        #[test]
        fn gcd_basic_examples() {
            // both even
            assert_eq!(
                gcd(&BigUint::from(48u32), &BigUint::from(18u32)),
                BigUint::from(6u32)
            );
            // both prime
            assert_eq!(
                gcd(&BigUint::from(17u32), &BigUint::from(13u32)),
                BigUint::from(1u32)
            );
            // equal numbers
            assert_eq!(
                gcd(&BigUint::from(16u32), &BigUint::from(16u32)),
                BigUint::from(16u32)
            );
            // one odd, one even
            assert_eq!(
                gcd(&BigUint::from(21u32), &BigUint::from(14u32)),
                BigUint::from(7u32)
            );
        }

        #[test]
        fn gcd_commutative() {
            let a = BigUint::from(123456u32);
            let b = BigUint::from(7890u32);

            assert_eq!(gcd(&a, &b), gcd(&b, &a));
        }

        #[test]
        fn gcd_matches_num_on_randoms() {
            use num_bigint::RandBigInt;
            use rand_chacha::{ChaCha20Rng, rand_core::SeedableRng};

            let mut rng = ChaCha20Rng::from_seed([42u8; 32]);

            for _ in 0..10 {
                let mut a = rng.gen_biguint(256);
                let mut b = rng.gen_biguint(256);

                // avoid case where both are zero
                if a.is_zero() && b.is_zero() {
                    a = BigUint::from(1u32);
                }

                let result = gcd(&a, &b);
                let reference = a.gcd(&b);
                assert_eq!(result, reference);
            }
        }
    }

    mod egcd_tests {
        use super::*;

        #[test]
        fn egcd_zero_cases() {
            let (gcd, _x, _y) = egcd(&BigInt::ZERO, &BigInt::ZERO);
            assert_eq!(gcd, BigInt::ZERO);

            let (gcd, _x, _y) = egcd(&BigInt::from(42), &BigInt::ZERO);
            assert_eq!(gcd, BigInt::from(42));

            let (gcd, _x, _y) = egcd(&BigInt::ZERO, &BigInt::from(85));
            assert_eq!(gcd, BigInt::from(85));
        }

        #[test]
        fn egcd_identity() {
            // a=240, b=46 -> gcd=2, and 240*x + 46*y = 2
            let a = BigInt::from(240);
            let b = BigInt::from(46);
            let (gcd, x, y) = egcd(&a, &b);

            // comparing to library gcd
            assert_eq!(gcd, a.gcd(&b));
            // check if bezout identity holds
            assert_eq!(gcd, a.clone() * x + b.clone() * y);
        }

        #[test]
        fn egcd_coprime() {
            let a = BigInt::from(17);
            let b = BigInt::from(13);
            let (gcd, x, y) = egcd(&a, &b);
            assert!(gcd.is_one());
            assert_eq!(gcd, a * x + b * y);
        }
    }

    mod modinv_tests {
        use super::*;

        #[test]
        fn modinv_basic_examples() {
            // 3^{-1} mod 11 = 4
            let inv = super::modinv(&BigInt::from(3), &BigInt::from(11)).unwrap();
            assert_eq!(inv, BigInt::from(4));

            // 10^{-1} mod 17 = 12
            let inv = super::modinv(&BigInt::from(10), &BigInt::from(17)).unwrap();
            assert_eq!(inv, BigInt::from(12));

            // 7^{-1} mod 26 = 15
            let inv = super::modinv(&BigInt::from(7), &BigInt::from(26)).unwrap();
            assert_eq!(inv, BigInt::from(15));
        }

        #[test]
        fn modinv_no_inverse_cases() {
            assert!(super::modinv(&BigInt::from(6), &BigInt::from(9)).is_none());
            assert!(super::modinv(&BigInt::from(12), &BigInt::from(18)).is_none());
        }

        #[test]
        fn modinv_handles_negative_a() {
            // (-3) mod 11 has inverse 7 (since -3 ≡ 8, and 8*7 ≡ 1 mod 11)
            let inv = super::modinv(&BigInt::from(-3), &BigInt::from(11)).unwrap();
            assert_eq!(inv, BigInt::from(7));
        }

        #[test]
        fn modinv_random_coprimes() {
            use num_bigint::RandBigInt;
            use num_traits::Signed;
            use rand_chacha::{ChaCha20Rng, rand_core::SeedableRng};

            // Small randomized check: pick random m, choose a coprime a,
            // verify (a * inv) % m == 1 with non-negative normalization.
            let mut rng = ChaCha20Rng::from_seed([42u8; 32]);

            for _ in 0..8 {
                // modulus in [2, 2^16)
                let m = rng.gen_bigint(16).abs() + BigInt::from(2);
                // choose a in [1, m-1]
                let mut a = rng.gen_bigint_range(&BigInt::one(), &m);

                // ensure coprime
                while a.gcd(&m) != BigInt::one() {
                    a = rng.gen_bigint_range(&BigInt::one(), &m);
                }

                let inv = super::modinv(&a, &m).expect("should be invertible");

                // Check (a * inv) % m == 1, normalized to [0, m)
                let one = (a * inv) % &m;
                let one = (one + &m) % &m;
                assert_eq!(one, BigInt::one());
            }
        }
    }
}
