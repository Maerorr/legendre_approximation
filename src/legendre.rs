use std::error::Error;

use crate::{Function, functions::{legendre_polynomial, horner, function_value}, integral::newton_cotes};

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

pub fn best_approximation(f: Function, eps: f64) -> (usize, bool) {
    let mut poly_deg = 1;
    let mut integral_nodes = 40;
    let mut results: Vec<(usize, f64)> = Vec::new();
    let mut error = eps + 1.;
    while error > eps {
        let lambdas = calculate_lambdas(f, poly_deg, integral_nodes);

        // part copied from main, that just calculates the error.
        let min = -1.;
        let max = 1.;
        let mut sum = 0.;
        let step = (max - min) / poly_deg as f64;
        for i in 0..poly_deg {
            sum += (function_value(min + i as f64 * step, f) - legendre_approx_value(&lambdas, min + i as f64 * step)).powi(2);
        }
        error = sum.sqrt();
        // end
        results.push((poly_deg, error));
        
        if poly_deg >= 10 {
            results.sort_by(
                |a, b| 
                a.1.partial_cmp(&b.1).unwrap());
            poly_deg = results[0].0;
            return (poly_deg, false);
        }

        poly_deg += 1;
        integral_nodes += 1;
    }
    (poly_deg, true)
}