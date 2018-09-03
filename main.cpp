/*
theta_t - a theta_xx + b kappa_a (theta^4 - phi) = q(t) f(x)
- alpha phi_xx + kappa_a (phi - theta^4) = 0,  x in (0, L), t in (0, T)
- a theta_x + beta (theta - theta_b) = 0 at x = 0
a theta_x + beta (theta - theta_b) = 0 at x = L
- alpha phi_x + gamma (phi - theta_b^4) = 0 at x = 0
alpha phi_x + gamma (phi - theta_b^4) = 0 at x = L
theta = theta_0 at t = 0

Problem: find q(t) for given r(t) = int_0^L g(x) theta(x,t) dx
*/

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <joker-fdm/bvp1d.h>

using namespace std;

// MODE_GIVEN_R - solve the inverse problem for given r(t)
// MODE_GIVEN_Q - calculate r(t) for given q(t), then solve the inverse problem
enum Mode {MODE_GIVEN_R, MODE_GIVEN_Q};
const Mode mode = MODE_GIVEN_Q;

const char* output_r_file_name = "output_r.txt";
const char* output_q_file_name = "output_q.txt";

// lengths of the space and time intervals
const double L = 50;
const double T = 30;
// coefficients in the equations and boundary conditions
const double a = 0.92;
const double alpha = 10.0/3;
const double kappa_a = 0.01;
const double b = 18.7;
const double beta = 10;
const double gamma = 0.3;
const double theta_b1 = 0.3;
const double theta_b2 = 0.8;

// initial condition theta_0(x)
double theta_0 (double x)
{
    const double theta_0_val = 1;
    return theta_0_val;
}

// function f(x) in the source term
double f (double x)
{
    if (x < L / 2)
        return 1.0e-1;
    else
        return 0.0;
}

// function g(x) in the integral
double g (double x)
{
    if (x > L / 2)
        return 1.0;
    else
        return 0.0;
}

// the given function q(t) in the source term for the case of mode == MODE_GIVEN_Q
double q_fun (double t)
{
    return (T - t) / T;
}

// the given function r(t) for the case of mode == MODE_GIVEN_R
double r_fun (double t)
{
    return 1.0;
}

// numbers of nodes of the grids
const int N = 50;
const int M = 50;

// parameters of the bisection method
const double q_init_guess = 0.0;  // initial guess
const double q_init_len = 0.1;  // initial length of the interval
const double q_tol = 1e-5;  // tolerance

double sign (double x)
{
    return (0 < x) - (x < 0);
}

// functions of the reaction terms and their derivatives

double theta_theta (double u)
{
    return b * kappa_a * pow(u, 4) * sign(u);
}

double d_theta_theta (double u)
{
    return b * kappa_a * 4 * pow(fabs(u), 3);
}

double theta_phi (double u)
{
    return - b * kappa_a * u;
}

double d_theta_phi (double u)
{
    return - b * kappa_a;
}

double phi_theta (double u)
{
    return - kappa_a * pow(u, 4) * sign(u);
}

double d_phi_theta (double u)
{
    return - kappa_a * 4 * pow(fabs(u), 3);
}

double phi_phi (double u)
{
    return kappa_a * u;
}

double d_phi_phi (double u)
{
    return kappa_a;
}

// calculate the solution at one time step
// sol is a vector of size 2:
//   initially sol contains the solution at the previous time step (m-1),
//   finally sol contains the solution at the current time step m
// q is the value of q(t) at the time interval (t_{m-1}, t_m)
// tau is the time grid step
void CalcSol (Data1D& data, vector<GridFunction1D>& sol, double q, double tau)
{
    int N = data.grid.K[0];
    for (int n = 0; n <= N; ++n) {
        double x = data.grid.coord(0, n);
        data.g[0](0, n) = q * f(x) + 1. / tau * sol[0](0, n);
        data.g[1](0, n) = 0.0;
    }
    SolveBVP1D(data, Parameters1D(), sol);
}

// calculate int_0^L g(x) theta(x) dx
double CalcIntegral (const Grid1D& grid, const GridFunction1D& theta)
{
    // trapezoid method
    double s = 0.0;
    for (int n = 0; n <= N; ++n) {
        double x = grid.coord(0, n);
        double theta_val = theta(0, n);
        double term = g(x) * theta_val;
        if (n == 0 || n == N)
            term /= 2;
        s += term;
    }
    return s * grid.h[0];
}

// sol1 and sol2 are vectors of size 2
// copy sol1 to sol2
void copy_sol(const Grid1D& grid,
    const vector<GridFunction1D>& sol1, vector<GridFunction1D>& sol2)
{
    for (int i = 0; i < 2; ++i) {
        for (int n = 0; n <= N; ++n)
            sol2[i](0, n) = sol1[i](0, n);
    }
}

int main ()
{
    vector<double> Lvec(1);  Lvec[0] = L;
    vector<int> Nvec(1);  Nvec[0] = N;
    Grid1D grid(Lvec, Nvec);
    Data1D data(2, grid);
    data.a[0][0] = a;  data.a[1][0] = alpha;
    data.b[0][0] = data.b[0][1] = beta;
    data.w[0][0] = beta * theta_b1;
    data.w[0][1] = beta * theta_b2;
    data.b[1][0] = data.b[1][1] = gamma;
    data.w[1][0] = gamma * pow(theta_b1, 4);
    data.w[1][1] = gamma * pow(theta_b2, 4);
    data.f[0][0][0] = theta_theta;  data.df[0][0][0] = d_theta_theta;
    data.f[0][0][1] = theta_phi;  data.df[0][0][1] = d_theta_phi;
    data.f[1][0][0] = phi_theta;  data.df[1][0][0] = d_phi_theta;
    data.f[1][0][1] = phi_phi;  data.df[1][0][1] = d_phi_phi;

    double tau = T / M;

    vector<double> c(2);
    c[0] = 1.0;  c[1] = 0.0;
    for (int i = 0; i < 2; ++i)
        data.c[i] = c[i] / tau;

    vector<GridFunction1D> sol(2), sol_prev(2);
    for (int i = 0; i < 2; ++i) {
        sol[i].set_grid(grid);
        sol_prev[i].set_grid(grid);
    }

    // r[m] = r(t_m) = int_0^L g(x) theta(x,t_m) dx, m = 0, 1, ..., M
    vector<double> r(M + 1);
    if (mode == MODE_GIVEN_R) {
        for (int m = 0; m <= M; ++m) {
            double t = m * tau;
            r[m] = r_fun(t);
        }
    }
    else {  // mode == MODE_GIVEN_Q
        // calculate r(t) for the given q(t)

        // set the function q(t)
        // q(t) = q[m], t in (t_{m-1}, t_m), m = 1, 2, ..., M
        vector<double> q(M + 1);
        for (int m = 0; m <= M; ++m) {
            double t = m * tau;
            q[m] = q_fun(t);
        }

        // set the initial condition
        {
            int i = 0;
            for (int n = 0; n <= N; ++n) {
                double x = grid.coord(0, n);
                sol[i](0, n) = theta_0(x);
            }
            // The solution for i = 1 (phi) is undefined.
            // The stationary equation for phi has to be solved
            //   in order to calculate phi at t = 0.
        }

        // solve the nonstationary problem and calculate r(t)
        // now sol[0] contains the initial function theta_0 and sol[1] is undefined
        r[0] = CalcIntegral(grid, sol[0]);
        for (int m = 1; m <= M; ++m) {
            // now sol contains the solution at the previous time step
            CalcSol(data, sol, q[m], tau);
            // now sol contains the solution at the current time step
            // calculate r(t_m) = int_0^L g(x) theta(x,t_m) dx
            r[m] = CalcIntegral(grid, sol[0]);
        }

        // output r(t)
        ofstream fout(output_r_file_name);
        fout.precision(10);
        for (int m = 0; m <= M; ++m)
            fout << tau * m << "   " << r[m] << endl;
    }  // if (mode)
    // r[m], m = 0, 1, ..., M, are calculated

    // solve the inverse problem - find q(t) for given r(t)

    // q(t) = q[m], t in (t_{m-1}, t_m), m = 1, 2, ..., M
    vector<double> q(M + 1);

    // sol_prev contains the solution at the previous time step
    // set the initial condition
    {
        int i = 0;
        for (int n = 0; n <= N; ++n) {
            double x = grid.coord(0, n);
            sol_prev[i](0, n) = theta_0(x);
        }
        // sol_prev[1] (phi) is undefined
    }

    double q_guess = q_init_guess;
    for (int m = 1; m <= M; ++m) {
        cout << "m = " << m << endl;
        // Denote by I(q) the value of the integral r(t_m) with q[m] = q.
        // Assume that I(q) is a monotonically increasing function.
        // find q[m] = q as the solution of the equation I(q) = r[m]

        // find q_1 and q_2 such that I(q_1) < r[m] and I(q_2) > r[m]
        double q_1, q_2, len1, len2, I;
        len1 = len2 = q_init_len / 4;  // length of the interval
        do {
            copy_sol(grid, sol_prev, sol);  // copy sol_prev to sol
            len1 *= 2;
            q_1 = q_guess - len1;
            CalcSol(data, sol, q_1, tau);
            I = CalcIntegral(grid, sol[0]);
        } while (I >= r[m]);
        do {
            copy_sol(grid, sol_prev, sol);  // copy sol_prev to sol
            len2 *= 2;
            q_2 = q_guess + len2;
            CalcSol(data, sol, q_2, tau);
            I = CalcIntegral(grid, sol[0]);
        } while (I <= r[m]);

        // apply the bisection method
        while (q_2 - q_1 >= q_tol) {
            q_guess = (q_1 + q_2) / 2;
            copy_sol(grid, sol_prev, sol);  // copy sol_prev to sol
            CalcSol(data, sol, q_guess, tau);
            I = CalcIntegral(grid, sol[0]);
            if (I < r[m])
                q_1 = q_guess;
            else
                q_2 = q_guess;
        }
        q[m] = q_guess;
        copy_sol(grid, sol, sol_prev);  // copy sol to sol_prev
    }

    // output q(t)
    ofstream fout(output_q_file_name);
    fout.precision(10);
    for (int m = 1; m <= M; ++m)
        fout << tau * m << "   " << q[m] << endl;

    return 0;
}
