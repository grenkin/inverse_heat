#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <joker-fdm/bvp1d.h>

using namespace std;

const char* output_r_file_name = "output_r.txt";

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
// initial condition
const double theta_init = 1;

// function f(x) in the source term
double f (double x)
{
    if (x < L / 2)
        return 1.0;
    else
        return 0.0;
}

// function in the integral
double g (double x)
{
    if (x > L / 2)
        return 1.0;
    else
        return 0.0;
}

// the given function q(t) in the source term
double q_fun (double t)
{
    return (T - t) / T;
}

// numbers of nodes of the grids
const int N = 50;
const int M = 50;

// functions of the reaction terms and their derivatives

double theta_theta (double u)
{
    return b * kappa_a * pow(u, 4);
}

double d_theta_theta (double u)
{
    return b * kappa_a * 4 * pow(u, 3);
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
    return - kappa_a * pow(u, 4);
}

double d_phi_theta (double u)
{
    return - kappa_a * 4 * pow(u, 3);
}

double phi_phi (double u)
{
    return kappa_a * u;
}

double d_phi_phi (double u)
{
    return kappa_a;
}

int main() {
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
    vector<double> c(2);
    c[0] = 1; c[1] = 0;
    double tau = T / M;
    for (int i = 0; i < 2; ++i)
        data.c[i] = c[i] / tau;

    vector<GridFunction1D> sol(2);
    vector<TimeGridFunction1D> sol_time(2);
    for (int i = 0; i < 2; ++i) {
        sol[i].set_grid(grid);
        sol_time[i].set_grid(grid, M);
    }

    // set the initial condition
    {
        int i = 0;
        for (int n = 0; n <= N; ++n) {
            sol[i](0, n) = theta_init;
            sol_time[i](0, 0, n) = sol[i](0, n);
        }
        // the solution for i = 1 (phi) is undefined
    }

    // set the function q(t)
    vector<double> q(M + 1);
    for (int m = 0; m <= M; ++m) {
        double t = m * tau;
        q[m] = q_fun(t);
    }

    // solve the nonstationary problem
    for (int m = 1; m <= M; ++m) {
        for (int n = 0; n <= N; ++n) {
            double x = grid.coord(0, n);
            data.g[0](0, n) = q[m] * f(x) + c[0] / tau * sol[0](0, n);
            data.g[1](0, n) = c[1] / tau * sol[1](0, n);
        }
        SolveBVP1D(data, Parameters1D(), sol);
        for (int i = 0; i < 2; ++i) {
            for (int n = 0; n <= N; ++n)
                sol_time[i](m, 0, n) = sol[i](0, n);
        }
    }

    // calculate the integral
    vector<double> r(M + 1);
    for (int m = 0; m <= M; ++m) {
        double s = 0.0;
        for (int n = 0; n <= N; ++n) {
            double x = grid.coord(0, n);
            double theta = sol_time[0](m, 0, n);
            double term = g(x) * theta;
            if (n == 0 || n == N)
                term /= 2;
            s += term;
        }
        r[m] = s * grid.h[0];
    }

    // output r(t)
    ofstream fout(output_r_file_name);
    for (int m = 0; m <= M; ++m)
        fout << tau * m << "   " << r[m] << endl;

    return 0;
}
