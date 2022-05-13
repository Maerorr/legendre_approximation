use crate::{Function, functions::{function_value, horner}};


// IMPLEMENTATION OF NEWTON COTES INTEGRATION FROM EXCERCISE 4.1

/// Calculates the numerator of the lambda coefficients
/// ### Integral of f(x) * L_k(x) from a to b
pub fn newton_cotes_top(f: Function, poly: Vec<f64>, a: f64, b: f64) -> f64 {
    let h = (b - a) / 2. as f64;
    let mut sum = 0.;
    sum += function_value(a, f) * horner(&poly, a);
    sum += 4. * function_value(a + h, f) * horner(&poly, a + h);
    sum += function_value(b, f) * horner(&poly, b);
    sum * h / 3.
}

/// Calculates the denominator of the lambda coefficients
/// ### Integral of L_k(x) * L_k(x) from a to b
pub fn newton_cotes_bot(_f: Function, poly: Vec<f64>, a: f64, b: f64) -> f64 {
    let h = (b - a) / 2. as f64;
    let mut sum = 0.;
    sum += horner(&poly, a) * horner(&poly, a);
    sum += 4. * horner(&poly, a+h) * horner(&poly, a+h);
    sum += horner(&poly, b) * horner(&poly, b);

    sum * h / 3.
}

/// Returns the value of the Newton-Cotes integration formula for the given function with given precision.
/// * f - chosen function from the Function enum
/// * poly - polynomial coefficient to be calculated using Horner's method
/// * nodes - amount of nodes for the Newton-Cotes integration
/// * which - changes between numerator and denominator of a lambda coefficient
pub fn newton_cotes(f: Function, poly: Vec<f64>, nodes: usize, which: bool) -> f64 
{
    let a = -1.;
    let b = 1.;

    match which {
        true => {
            let h = (b - a) / (nodes as f64);
            let mut sum = 0.;
            let mut x = a;
            for i in 0..nodes {
                sum += newton_cotes_top(f, poly.to_vec(), x, x + h);
                x += h;
            }
            return sum
        },
        false => {
            let h = (b - a) / (nodes as f64);
            let mut sum = 0.;
            let mut x = a;
            for i in 0..nodes {
                sum += newton_cotes_bot(f, poly.to_vec(), x, x + h);
                x += h;
            }
            return sum
        },
    };  
}