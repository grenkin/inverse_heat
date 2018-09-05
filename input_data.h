#ifndef INPUT_DATA_H_INCLUDED
#define INPUT_DATA_H_INCLUDED

#include <vector>

// MODE_GIVEN_R - solve the inverse problem for given r(t)
// MODE_GIVEN_Q - calculate r(t) for given q(t), then solve the inverse problem
enum Mode {MODE_GIVEN_R, MODE_GIVEN_Q};

struct InputData {
    // lengths of the space and time intervals
    double L, T;
    // coefficients in the equations and boundary conditions
    double a, alpha, kappa_a, b, beta, gamma, theta_b1, theta_b2;
    // mode: whether r(t) is given or q(t) is given
    Mode mode;
    // initial condition theta_0(x)
    std::vector<double> theta_0;  // vector of length (N + 1)
    // function f(x) in the source term
    std::vector<double> f;  // vector of length (N + 1)
    // function g(x) in the integral
    std::vector<double> g;  // vector of length (N + 1)
    // the given function q(t) in the source term
    // if mode == MODE_GIVEN_Q then q(t) is optional
    std::vector<double> q;  // vector of length (M + 1)
    // the given function r(t)
    // if mode == MODE_GIVEN_R then r(t) is optional
    std::vector<double> r;  // vector of length (M + 1)
    // parameters of the bisection method
    double q_init_guess;  // initial guess
    double q_init_len;  // initial length of the interval
    double q_tol;  // tolerance
    // parameters of verification of monotonicity of the integral depending on q
    bool verify_monotonicity;
    // q is set from q_1 to q_2 with step q_step
    // if verify_monotonicity == false then the following parameters are optional
    double monot_ver_q_1, monot_ver_q_2, monot_ver_q_step;
    // numbers of nodes of the grids: n = 0, ..., N - space, m = 0, ..., M - time
    int N, M;

    InputData (std::string);
};

#endif // INPUT_DATA_H_INCLUDED
