#include "input_data.h"

double theta_0_fun (double x)
{
    const double theta_0_val = 1;
    return theta_0_val;
}

double f_fun (double x, double L)
{
    if (x < L / 2)
        return 1.0e-1;
    else
        return 0.0;
}

double g_fun (double x, double L)
{
    if (x > L / 2)
        return 1.0;
    else
        return 0.0;
}

double q_fun (double t, double T)
{
    return (T - t) / T;
}

double r_fun (double t, double T)
{
    return 20.0 + 5 * (T - t) / T;
}

InputData::InputData ()
{
    L = 50;
    T = 30;
    a = 0.92;
    alpha = 10.0/3;
    kappa_a = 0.01;
    b = 18.7;
    beta = 10;
    gamma = 0.3;
    theta_b1 = 0.3;
    theta_b2 = 0.8;
    mode = MODE_GIVEN_Q;
    q_init_guess = 0.0;
    q_init_len = 0.08;
    q_tol = 1e-5;
    verify_monotonicity = true;
    monot_ver_q_1 = -10;
    monot_ver_q_2 = 10 + 1e-5;
    monot_ver_q_step = 0.1;
    N = 50;
    M = 50;
    double h = L / N, tau = T / M;
    theta_0.resize(N + 1);
    f.resize(N + 1);
    g.resize(N + 1);
    for (int n = 0; n <= N; ++n) {
        double x = h * n;
        theta_0[n] = theta_0_fun(x);
        f[n] = f_fun(x, L);
        g[n] = g_fun(x, L);
    }
    q.resize(M + 1);
    r.resize(M + 1);
    for (int m = 0; m <= M; ++m) {
        double t = tau * m;
        q[m] = q_fun(t, T);
        r[m] = r_fun(t, T);
    }
}
