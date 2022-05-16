use eframe::{
    egui::{self, plot::{Points, Plot, Values, Value, Line, VLine}, Layout},
    epi::{App}, run_native,
};
use functions::function_value;
use legendre::*;

mod functions;
mod legendre;
mod integral;

#[derive(Clone, Copy, PartialEq)]
pub enum Mode {
    Nodes,
    AproxError,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Function {
    Poly1,
    Poly2,
    PerfectFit,
    Linear,
    Sinusoidal,
    Absolute,
    Mixed,
}

struct AppState {
    function: Function,
    no_of_nodes: usize,
    mode: Mode,
    chosen_function_values: Vec<Value>,
    approx_values: Vec<Value>,
    lambdas: Vec<f64>,
    center_plot: bool,
    integral_nodes: usize,
    approx_error: f64,
    polynomial: String,
    given_approx_error: f64,
    epsilon_flag: bool,
}

impl AppState {
    fn new() -> AppState {
        AppState {
            function: Function::Poly1,
            no_of_nodes: 2,
            mode: Mode::Nodes,
            chosen_function_values: Vec::new(),
            approx_values: Vec::new(),
            lambdas: Vec::new(),
            center_plot: false,
            integral_nodes: 2,
            approx_error: 0.,
            polynomial: String::new(),
            given_approx_error: 0.1,
            epsilon_flag: true,
        }
    }
}

impl App for AppState {
    fn name(&self) -> &str {
        "Laguere Polynomial Approximation"
    }

    fn update(&mut self, ctx: &eframe::egui::Context, _frame: &eframe::epi::Frame) {
        ctx.set_pixels_per_point(1.5);
        egui::SidePanel::left("left_panel").min_width(150.).show(ctx, |ui| {
            egui::ScrollArea::vertical().show(ui, |ui| {

                // ##################################
                //         FUNCTION SELECTION
                // ##################################

                ui.group(|ui| {
                    ui.heading("Function");
                    ui.add_space(5.);
    
                    ui.with_layout(Layout::top_down(egui::Align::LEFT), |ui| {
                        ui.radio_value(&mut self.function, Function::Poly1, "Polynomial 2nd Power");
                        ui.radio_value(&mut self.function, Function::Poly2, "Polynomial 4th power");
                        ui.radio_value(&mut self.function, Function::PerfectFit, "Perfect Fit");
                        ui.radio_value(&mut self.function, Function::Linear, "Linear");
                        ui.radio_value(&mut self.function, Function::Sinusoidal, "Sinusoidal");
                        ui.radio_value(&mut self.function, Function::Absolute, "Absolute");
                        ui.radio_value(&mut self.function, Function::Mixed, "Mixed");
                    });
                });

                // ##################################
                //         APPROX. RANGE
                // ##################################
                if ui.button("Nodes").clicked() {
                    self.mode = Mode::Nodes;
                }
                if ui.button("Approx. Error").clicked() {
                    self.mode = Mode::AproxError;
                }
                match self.mode {
                    Mode::Nodes => {
                        ui.group(|ui| {
                            ui.group(|ui| {
                                ui.label("Approximation takes place on the interval [-1, 1]");
                            }); 
                            //ui.label("Mode");
                            //ui.radio_value(&mut self.mode, Mode::Nodes, "Nodes");
                            //ui.radio_value(&mut self.mode, Mode::AproxError, "Approx. Error");
                            ui.group(|ui| {
                                ui.label("Polynomial Degree");
                                ui.add(egui::Slider::new(&mut self.no_of_nodes, 2..=10));
                                ui.label("Newton-Cotes Nodes");
                                ui.add(egui::Slider::new(&mut self.integral_nodes, 2..=40));
                            });
                            if ui.button("Calculate").clicked() {
                                if self.integral_nodes < self.no_of_nodes {
                                    self.integral_nodes = self.no_of_nodes;
                                }
                                // default some parameters
                                self.chosen_function_values = Vec::new();
                                self.approx_values = Vec::new();
                                self.lambdas = Vec::new();
            
                                let min = -1.;
                                let max = 1.;
                                // generating values of chosen function for the plot
                                self.chosen_function_values = (0..10000)
                                .map(|i| {
                                    let x = min + (i as f64 *
                                    ((max) - (min)) / 10000.);
                                    Value::new(x, function_value(x, self.function))
                                })
                                .collect();
            
                                // generating values of approximated function for the plot
                                self.lambdas = calculate_lambdas(self.function, self.no_of_nodes, self.integral_nodes);
                                self.approx_values = (0..10000)
                                .map(|i| {
                                    let x = min + (i as f64 *
                                    ((max) - (min)) / 10000.);
                                    Value::new(x, legendre_approx_value(&self.lambdas, x))
                                })
                                .collect();
            
                                let mut sum = 0.;
                                let step = (max - min) / self.no_of_nodes as f64;
                                for i in 0..self.no_of_nodes {
                                    sum += (function_value(min + i as f64 * step, self.function) - legendre_approx_value(&self.lambdas, min + i as f64 * step)).powi(2);
                                }
                                self.approx_error = sum.sqrt();
                                
                                let mut polynomial: String = String::from(" ");
                                let poly = get_coefficients(&self.lambdas);
                                for (i, j) in poly.iter().enumerate() {
                                    if i == self.no_of_nodes  {
                                        polynomial += format!("{:.3}x^{}" , j, self.no_of_nodes - i).as_str();
                                    } else {
                                        polynomial += format!("{:.3}x^{} + " , j, self.no_of_nodes - i).as_str();
                                    }
                                    
                                }
                                self.polynomial = polynomial;
                            }
                            let error = format!("Approx. Error: {:.6}", self.approx_error);
                            ui.group(|ui| {
                                ui.label(error);
                            });
                        });
                    },
                    Mode::AproxError => {
                        ui.group(|ui| {
                            self.integral_nodes = 40;
                            ui.group(|ui| {
                                ui.label("Approximation takes place on the interval [-1, 1]");
                            }); 
                            ui.group(|ui| {
                                ui.label("Approx. Epsilon: ");
                                ui.add(egui::Slider::new(&mut self.given_approx_error, 1e-15..=0.1).logarithmic(true));
                            });
                            if ui.button("Calculate").clicked() {
                                // default some parameters
                                self.chosen_function_values = Vec::new();
                                self.approx_values = Vec::new();
                                self.lambdas = Vec::new();
            
                                let (best_deg, flag) = best_approximation(self.function, self.given_approx_error);
                                self.no_of_nodes = best_deg;
                                self.epsilon_flag = flag;

                                let min = -1.;
                                let max = 1.;
                                // generating values of chosen function for the plot
                                self.chosen_function_values = (0..10000)
                                .map(|i| {
                                    let x = min + (i as f64 *
                                    ((max) - (min)) / 10000.);
                                    Value::new(x, function_value(x, self.function))
                                })
                                .collect();
            
                                // generating values of approximated function for the plot
                                self.lambdas = calculate_lambdas(self.function, self.no_of_nodes, self.integral_nodes);
                                self.approx_values = (0..10000)
                                .map(|i| {
                                    let x = min + (i as f64 *
                                    ((max) - (min)) / 10000.);
                                    Value::new(x, legendre_approx_value(&self.lambdas, x))
                                })
                                .collect();
            
                                let mut sum = 0.;
                                let step = (max - min) / self.no_of_nodes as f64;
                                for i in 0..self.no_of_nodes {
                                    sum += (function_value(min + i as f64 * step, self.function) - legendre_approx_value(&self.lambdas, min + i as f64 * step)).powi(2);
                                }
                                self.approx_error = sum.sqrt();
                                
                                let mut polynomial: String = String::from(" ");
                                let poly = get_coefficients(&self.lambdas);
                                for (i, j) in poly.iter().enumerate() {
                                    if i == self.no_of_nodes  {
                                        polynomial += format!("{:.3}x^{}" , j, self.no_of_nodes - i).as_str();
                                    } else {
                                        polynomial += format!("{:.3}x^{} + " , j, self.no_of_nodes - i).as_str();
                                    }
                                    
                                }
                                self.polynomial = polynomial;
                            }
                            let error = format!("Approx. Error: {:.6}", self.approx_error);
                            ui.group(|ui| {
                                if !self.epsilon_flag {
                                    ui.label("The search criteria were not met. Displaying the best approximation found.");
                                }
                                ui.label(error);
                            });
                        });
                    },
                }
                
                
            });
        });
        egui::CentralPanel::default().show(ctx, |ui| {
            egui::ScrollArea::vertical().show(ui, |ui| {
                //chosen function
                let chosen_values = Values::from_values(self.chosen_function_values.clone());
                let approximated_values = Values::from_values(self.approx_values.clone());
                let chosen_plot = Line::new(chosen_values).name("Chosen Function");
                let approx_plot = Line::new(approximated_values).name("Approx. Function");

                let vline_left = VLine::new(-1.);
                let vline_right = VLine::new(1.);

                ui.checkbox(&mut self.center_plot, "Center Plot");
                let mut plot = Plot::new("my_plot")
                    .show_x(true)
                    .show_y(true)
                    .legend(egui::widgets::plot::Legend::default());
                plot = plot.data_aspect(1.0);
                if self.center_plot {
                    plot = plot.center_x_axis(true).center_y_axis(true)
                }
                plot.show(ui, |plot_ui| {
                    plot_ui.line(chosen_plot);
                    plot_ui.line(approx_plot);
                    plot_ui.vline(vline_left);
                    plot_ui.vline(vline_right);
                });
                ui.group(|ui| {
                    ui.add_space(5.);
                    ui.label(self.polynomial.as_str());
                    ui.add_space(5.);
                });
            }); 
        });
    }
}

fn main() {
    let native_options = eframe::NativeOptions {
        initial_window_size: Some(egui::vec2(1000., 660.0)),
        ..eframe::NativeOptions::default()
    };
    run_native(Box::new(AppState::new()), native_options);
}
