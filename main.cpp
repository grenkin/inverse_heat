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
#include <cstdlib>
#include <conio.h>
#include <joker-fdm/bvp1d.h>
#include "input_data.h"

using namespace std;

const char* input_file_name = "input.txt";
const char* output_r_file_name = "output_r.txt";
const char* output_q_file_name = "output_q.txt";
const char* output_log_file_name = "output_log.txt";
const char* output_monot_log_file_name = "output_monot_log.txt";

// lengths of the space and time intervals
double L;
double T;
// coefficients in the equations and boundary conditions
double a;
double alpha;
double kappa_a;
double b;
double beta;
double gamma;
double theta_b1;
double theta_b2;
// numbers of nodes of the grids
int N;
int M;

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
void CalcSol (Data1D& data, vector<GridFunction1D>& sol, double q,
    const InputData& id)
{
    int N = data.grid.K[0];
    double tau = id.T / id.M;
    for (int n = 0; n <= N; ++n) {
        data.g[0](0, n) = q * id.f[n] + 1. / tau * sol[0](0, n);
        data.g[1](0, n) = 0.0;
    }
    Parameters1D param;
    param.sol_method = SOL_METHOD_UMFPACK;
    param.Newton_tol = 1e-10;
    SolveBVP1D(data, param, sol);
}

// calculate int_0^L g(x) theta(x) dx
double CalcIntegral (const Grid1D& grid, const GridFunction1D& theta,
    const InputData& id)
{
    // trapezoid method
    double s = 0.0;
    for (int n = 0; n <= N; ++n) {
        double theta_val = theta(0, n);
        double term = id.g[n] * theta_val;
        if (n == 0 || n == N)
            term /= 2;
        s += term;
    }
    return s * grid.h[0];
}

// sol1 and sol2 are vectors of size 2
// copy sol1 to sol2
void copy_sol (const Grid1D& grid,
    const vector<GridFunction1D>& sol1, vector<GridFunction1D>& sol2)
{
    for (int i = 0; i < 2; ++i) {
        for (int n = 0; n <= N; ++n)
            sol2[i](0, n) = sol1[i](0, n);
    }
}

int main ()
{
    InputData id(input_file_name);
    L = id.L;
    T = id.T;
    a = id.a;
    alpha = id.alpha;
    kappa_a = id.kappa_a;
    b = id.b;
    beta = id.beta;
    gamma = id.gamma;
    theta_b1 = id.theta_b1;
    theta_b2 = id.theta_b2;
    N = id.N;
    M = id.M;

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
    if (id.mode == MODE_GIVEN_R) {
        for (int m = 0; m <= M; ++m)
            r[m] = id.r[m];
    }
    else {  // mode == MODE_GIVEN_Q
        // calculate r(t) for the given q(t)

        cout << "Calculate r(t)..." << endl;

        // set the function q(t)
        // q(t) = q[m], t in (t_{m-1}, t_m), m = 1, 2, ..., M
        vector<double> q(M + 1);
        for (int m = 1; m <= M; ++m)
            q[m] = id.q[m];

        // set the initial condition
        {
            int i = 0;
            for (int n = 0; n <= N; ++n)
                sol[i](0, n) = id.theta_0[n];
            // The solution for i = 1 (phi) is undefined.
            // The stationary equation for phi has to be solved
            //   in order to calculate phi at t = 0.
        }

        // solve the nonstationary problem and calculate r(t)
        // now sol[0] contains the initial function theta_0 and sol[1] is undefined
        r[0] = CalcIntegral(grid, sol[0], id);
        for (int m = 1; m <= M; ++m) {
            // now sol contains the solution at the previous time step
            CalcSol(data, sol, q[m], id);
            // now sol contains the solution at the current time step
            // calculate r(t_m) = int_0^L g(x) theta(x,t_m) dx
            r[m] = CalcIntegral(grid, sol[0], id);
        }

        // output r(t)
        ofstream fout(output_r_file_name);
        fout.precision(10);
        for (int m = 0; m <= M; ++m)
            fout << tau * m << "   " << r[m] << endl;
    }  // if (mode)
    // r[m], m = 0, 1, ..., M, are calculated

    // solve the inverse problem - find q(t) for given r(t)

    ofstream flog(output_log_file_name);
    ofstream flog_monot(output_monot_log_file_name);
    flog_monot.precision(10);
    flog_monot << "q" << endl;
    for (double q = id.monot_ver_q_1; q <= id.monot_ver_q_2;
        q += id.monot_ver_q_step)
    {
        flog_monot << q << "  ";
    }
    flog_monot << "\n\n\nI(q)\n\n";

    // q(t) = q[m], t in (t_{m-1}, t_m), m = 1, 2, ..., M
    vector<double> q(M + 1);

    // sol_prev contains the solution at the previous time step
    // set the initial condition
    {
        int i = 0;
        for (int n = 0; n <= N; ++n)
            sol_prev[i](0, n) = id.theta_0[n];
        // sol_prev[1] (phi) is undefined
    }

    double q_guess = id.q_init_guess;
    double q_len = id.q_init_len / 4;
    for (int m = 1; m <= M; ++m) {
        cout << "m = " << m << endl;
        // Denote by I(q) the value of the integral r(t_m) with q[m] = q.
        // Assume that I(q) is a monotonically increasing function.

        if (id.verify_monotonicity) {
            // verify monotonicity of the function I(q)
            cout << "Verify monotonicity... ";
            bool monotone = true;
            flog_monot << "m = " << m << endl;
            bool start = true;
            double I_last;
            for (double q = id.monot_ver_q_1; q <= id.monot_ver_q_2;
                q += id.monot_ver_q_step)
            {
                copy_sol(grid, sol_prev, sol);  // copy sol_prev to sol
                CalcSol(data, sol, q, id);
                double I = CalcIntegral(grid, sol[0], id);
                if (monotone && !start && I_last > I) {
                    cout << "Monotonicity of I(q) is not fulfilled!!!\n";
                    cout << "q = " << q << endl;
                    monotone = false;
                }
                start = false;
                I_last = I;
                flog_monot << I << "  ";
            }
            if (!monotone) {
                // getch();
                // flog_monot.close();
                // exit(1);
            }
            else
                cout << "OK" << endl;
            flog_monot << "\n";
        }

        // find q[m] = q as the solution of the equation I(q) = r[m]

        flog << "m = " << m << endl;
        // now q_guess contains q(t) from the previous time step
        double q_guess_old = q_guess;
        // find q_1 and q_2 such that I(q_1) < r[m] and I(q_2) > r[m]
        double q_1, q_2, len1, len2, I;
        len1 = len2 = q_len;  // length of the interval
        flog << "len1 =";
        do {
            copy_sol(grid, sol_prev, sol);  // copy sol_prev to sol
            len1 *= 2;
            q_1 = q_guess - len1;
            CalcSol(data, sol, q_1, id);
            I = CalcIntegral(grid, sol[0], id);
            flog << "  " << len1;
        } while (I >= r[m]);
        flog << "\nlen2 =";
        do {
            copy_sol(grid, sol_prev, sol);  // copy sol_prev to sol
            len2 *= 2;
            q_2 = q_guess + len2;
            CalcSol(data, sol, q_2, id);
            I = CalcIntegral(grid, sol[0], id);
            flog << "  " << len2;
        } while (I <= r[m]);
        flog << "\nq_1 = " << q_1 << "   q_2 = " << q_2 << endl;

        // apply the bisection method
        flog << "q_guess =";
        while (q_2 - q_1 >= id.q_tol) {
            q_guess = (q_1 + q_2) / 2;
            copy_sol(grid, sol_prev, sol);  // copy sol_prev to sol
            CalcSol(data, sol, q_guess, id);
            I = CalcIntegral(grid, sol[0], id);
            if (I < r[m])
                q_1 = q_guess;
            else
                q_2 = q_guess;
            flog << "  " << q_guess;
        }
        q[m] = q_guess;
        copy_sol(grid, sol, sol_prev);  // copy sol to sol_prev
        flog << "\n\n";
        // now q_guess contains q(t) at the current time step
        q_len = fabs(q_guess - q_guess_old);
    }

    // output q(t)
    ofstream fout(output_q_file_name);
    fout.precision(10);
    for (int m = 1; m <= M; ++m)
        fout << tau * m << "   " << q[m] << endl;

    cout << "Done";
    getch();
    return 0;
}
