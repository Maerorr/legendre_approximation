use std::{
    f64::consts::{E, PI},
};

use crate::Function;

pub fn polynomial1(x: f64) -> f64 {
    // 0.15x^2 - x - 1
    -1. + x * (-1. + x * 0.15)
}

pub fn polynomial2(x: f64) -> f64 {
    //0.07*x^4+0.3*x^3-0.2*x^2-x-1.
    -1. + x * (-1. + x * (-0.2 + x * (-0.3 + x * 0.07)))
    //2.*E.powf(-3.*x)*(4.*x + 3.*x*x)
}

pub fn perfect_fit(x: f64) -> f64 {
    PI * (x - 3.).powi(-2) + 0.1
}

pub fn linear(x: f64) -> f64 {
    0.5 * x + 2.
}

pub fn sinusoidal(x: f64) -> f64 {
    x.sin()
}

pub fn absolute(x: f64) -> f64 {
    //((x-2.).abs()-2.).abs()
    x.abs()
}

pub fn mixed(x: f64) -> f64 {
    ((x - 2.).abs() - 2.).abs() + x.sin() + 0.05 * x.powf(3.)
}

pub fn function_value(x: f64, func: Function) -> f64 {
    match func {
        Function::Poly1 => polynomial1(x),
        Function::Poly2 => polynomial2(x),
        Function::PerfectFit => perfect_fit(x),
        Function::Linear => linear(x),
        Function::Sinusoidal => sinusoidal(x),
        Function::Absolute => absolute(x),
        Function::Mixed => mixed(x),
    }
}

pub fn factorial(n: usize) -> i64 {
    if n == 0 {
        1
    } else {
        n as i64 * factorial(n - 1)
    }
}

pub fn binomial_coeff(top: usize, bot: usize) -> i64 {
    if top > bot {
        factorial(top) / (factorial(top - bot) * factorial(bot))
    }
    else if top == bot {
        1
    } else {
        0
    }
}

pub fn pow(x: i64, n: usize) -> i64 {
    if n == 0 {
        return 1
    }
    if n == 1 {
        return x
    }
    let mut out = x;
    for i in 1..n {
        out *= x;
    }
    out
}

/// Returns the value of a function in point x.
/// Uses Horner's method.
/// * a - vector of coefficients of a function for example 3x^2 + 2x + 1 = {1, 2, 3}
pub fn horner(a: &Vec<f64>, x: f64) -> f64 {
    let mut i = 0;
    let mut out = a[i as usize];
    i += 1;
    while i < a.len() {
        out *= x;
        out += a[i as usize];
        i += 1;
    }
    out
}

/// https://en.wikipedia.org/wiki/Legendre_polynomials#Rodrigues'_formula_and_other_explicit_formulas
pub fn legendre_polynomial(deg: usize) -> Vec<f64> {
    let mut out: Vec<f64> = Vec::new();

    for k in 0..(deg/2)+1 {
        out.push( ((-1i64).pow(k as u32) * binomial_coeff(deg, k) * binomial_coeff(2*deg - 2*k, deg)) as f64);
        out.push(0.);
    }
    if deg % 2 != 1 {
        out.pop();
    }
    
    out
}