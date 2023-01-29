// -*- mode: rust; -*-
//
// This file is part of curve25519-dalek.
// Copyright (c) 2016-2021 isis lovecruft
// Copyright (c) 2016-2019 Henry de Valence
// See LICENSE for licensing information.
//
// Authors:
// - isis agora lovecruft <isis@patternsinthevoid.net>
// - Henry de Valence <hdevalence@hdevalence.ca>
#![allow(non_snake_case)]

use core::time::Duration;

use backend::serial::curve_models::{ProjectiveNielsPoint, ProjectivePoint};
use constants;
use edwards::EdwardsPoint;
use scalar::Scalar;
use traits::Identity;
use window::NafLookupTable5;

/// Compute \\(aA + bB\\) in variable time, where \\(B\\) is the Ed25519 basepoint.
pub fn mul(a: &Scalar, A: &EdwardsPoint, b: &Scalar) -> EdwardsPoint {
    let a_naf = a.non_adjacent_form(5);
    let b_naf = b.non_adjacent_form(8);

    // Find starting index
    let mut i: usize = 255;
    for j in (0..256).rev() {
        i = j;
        if a_naf[i] != 0 || b_naf[i] != 0 {
            break;
        }
    }

    let table_A = NafLookupTable5::<ProjectiveNielsPoint>::from(A);
    let table_B = &constants::AFFINE_ODD_MULTIPLES_OF_BASEPOINT;

    let mut r = ProjectivePoint::identity();
    loop {
        let mut t = r.double();

        if a_naf[i] > 0 {
            t = &t.to_extended() + &table_A.select(a_naf[i] as usize);
        } else if a_naf[i] < 0 {
            t = &t.to_extended() - &table_A.select(-a_naf[i] as usize);
        }

        if b_naf[i] > 0 {
            t = &t.to_extended() + &table_B.select(b_naf[i] as usize);
        } else if b_naf[i] < 0 {
            t = &t.to_extended() - &table_B.select(-b_naf[i] as usize);
        }

        r = t.to_projective();

        if i == 0 {
            break;
        }
        i -= 1;
    }

    r.to_extended()
}

/// Compute \\(aA + bB\\) in variable time, where \\(B\\) is the Ed25519 basepoint.
pub fn mul_timed(
    a: &Scalar,
    A: &EdwardsPoint,
    b: &Scalar,
) -> (EdwardsPoint, Duration, Duration, Duration) {
    let a_naf = a.non_adjacent_form(5);
    let b_naf = b.non_adjacent_form(8);
    println!("a_naf: {a_naf:?}");
    println!("b_naf: {b_naf:?}");

    // Find starting index
    let mut i: usize = 255;
    for j in (0..256).rev() {
        i = j;
        if a_naf[i] != 0 || b_naf[i] != 0 {
            break;
        }
    }
    let mut a = std::time::Duration::ZERO;
    let mut b = std::time::Duration::ZERO;

    let clock = std::time::Instant::now();
    let table_A = NafLookupTable5::<ProjectiveNielsPoint>::from(A);
    let table_B = &constants::AFFINE_ODD_MULTIPLES_OF_BASEPOINT;
    let table = clock.elapsed();

    let mut r = ProjectivePoint::identity();
    let mut byzantine_counter = 0;
    let mut byzantine_inefficiency_counter = 0;
    loop {
        byzantine_counter += 1;
        let mut t = r.double();

        let clock = std::time::Instant::now();
        if a_naf[i] > 0 {
            t = &t.to_extended() + &table_A.select(a_naf[i] as usize);
        } else if a_naf[i] < 0 {
            t = &t.to_extended() - &table_A.select(-a_naf[i] as usize);
        } else {
            // println!("a={i}");
            byzantine_inefficiency_counter += 1;
        }
        a += clock.elapsed();

        let clock = std::time::Instant::now();
        if b_naf[i] > 0 {
            t = &t.to_extended() + &table_B.select(b_naf[i] as usize);
        } else if b_naf[i] < 0 {
            t = &t.to_extended() - &table_B.select(-b_naf[i] as usize);
        } else {
            // println!("b={i}");
            byzantine_inefficiency_counter += 1;
        }
        b += clock.elapsed();

        r = t.to_projective();

        if i == 0 {
            break;
        }
        i -= 1;
    }
    println!("iters={byzantine_counter} inefficiency={byzantine_inefficiency_counter} byzantine_score={}", 2*byzantine_counter-byzantine_inefficiency_counter);

    (r.to_extended(), table, a, b)
}

/// Compute \\(aA + bB\\) in variable time, where \\(B\\) is the Ed25519 basepoint.
pub fn mul_byz_score(
    a: &Scalar,
    A: &EdwardsPoint,
    b: &Scalar,
) -> (EdwardsPoint, usize) {
    let a_naf = a.non_adjacent_form(5);
    let b_naf = b.non_adjacent_form(8);

    // Find starting index
    let mut i: usize = 255;
    for j in (0..256).rev() {
        i = j;
        if a_naf[i] != 0 || b_naf[i] != 0 {
            break;
        }
    }

    let table_A = NafLookupTable5::<ProjectiveNielsPoint>::from(A);
    let table_B = &constants::AFFINE_ODD_MULTIPLES_OF_BASEPOINT;

    let mut r = ProjectivePoint::identity();
    let mut byzantine_counter = 0;
    let mut byzantine_inefficiency_counter = 0;
    loop {
        byzantine_counter += 1;
        let mut t = r.double();

        if a_naf[i] > 0 {
            t = &t.to_extended() + &table_A.select(a_naf[i] as usize);
        } else if a_naf[i] < 0 {
            t = &t.to_extended() - &table_A.select(-a_naf[i] as usize);
        } else {
            byzantine_inefficiency_counter += 1;
        }

        if b_naf[i] > 0 {
            t = &t.to_extended() + &table_B.select(b_naf[i] as usize);
        } else if b_naf[i] < 0 {
            t = &t.to_extended() - &table_B.select(-b_naf[i] as usize);
        } else {
            byzantine_inefficiency_counter += 1;
        }

        r = t.to_projective();

        if i == 0 {
            break;
        }
        i -= 1;
    }
    let score = 2*byzantine_counter-byzantine_inefficiency_counter;

    (r.to_extended(), score)
}
