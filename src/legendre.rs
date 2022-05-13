use crate::{Function, functions::{legendre_polynomial, horner}, integral::newton_cotes};

/// calculates lambdas for the approximation polynomial
/// * f - Function from the Function enum
/// * poly_deg - Degree of the approximating polynomial
/// * integral_nodes - Amount of nodes for the Newton-Cotes integral
pub fn calculate_lambdas(f: Function, poly_deg: usize, integral_nodes: usize) -> Vec<f64> {
    let mut out: Vec<f64> = Vec::new();
    for i in 0..(poly_deg+1) {
        let poly = legendre_polynomial(i);
        out.push(
            newton_cotes(f, poly.to_vec(), integral_nodes , true) /
            newton_cotes(f, poly.to_vec(), integral_nodes , false) 
        )
    }
    out
}

pub fn legendre_approx_value(lambdas: &Vec<f64> , x: f64) -> f64 {
    let mut sum = 0.;
    for (i, lambda) in lambdas.iter().enumerate() {
        sum += lambda * horner(&legendre_polynomial(i), x);
    }
    sum
}

pub fn get_coefficients(lambdas: &Vec<f64>) -> Vec<f64> {
    let mut out: Vec<f64> = vec![0.; lambdas.len()];
    for i in 0..lambdas.len() {
        let mut temp: Vec<f64> = Vec::new();
        let poly = legendre_polynomial(i);
        for j in 0..poly.len() {
            temp.push(poly[j] * lambdas[i]);
        }
        for (iter, elem) in temp.iter().rev().enumerate() {
            out[lambdas.len() - 1 - iter] += elem;
        }
    }
    out
}