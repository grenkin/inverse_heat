#include <iostream>
#include <cstdio>
#include <fstream>
#include <conio.h>
#include <boost/program_options.hpp>
#include "input_data.h"

namespace po = boost::program_options;
using namespace std;

struct ParamIsNotSet {
    std::string param;
    ParamIsNotSet (std::string param)
        : param(param)
    {}
};

struct ParamIsNotPositive {
    std::string param;
    ParamIsNotPositive (std::string param)
        : param(param)
    {}
};

void print_error (const string msg)
{
    cerr << msg << "\n";
    getch();
    exit(1);
}

void get_double_param (const po::variables_map& vm, const char* param,
    double& val)
{
    if (vm.count(param))
        val = vm[param].as<double>();
    else
        throw ParamIsNotSet(param);
}

void get_pos_double_param (po::variables_map& vm, const char* param,
    double& val)
{
    if (vm.count(param)) {
        val = vm[param].as<double>();
        if (val < 0)
            throw ParamIsNotPositive(param);
    }
    else
        throw ParamIsNotSet(param);
}

void get_int_param (const po::variables_map& vm, const char* param, int& val)
{
    if (vm.count(param)) {
        val = vm[param].as<int>();
        if (val < 0)
            throw ParamIsNotPositive(param);
    }
    else
        throw ParamIsNotSet(param);
}

void get_string_param (const po::variables_map& vm, const char* param,
    string& val)
{
    if (vm.count(param))
        val = vm[param].as<string>();
    else
        throw ParamIsNotSet(param);
}

template <typename T>
std::string toString (T val)
{
    std::ostringstream oss;
    oss << val;
    return oss.str();
}

extern "C" FILE* popen(const char* command, const char* mode);
extern "C" int pclose(FILE* stream);

// call program cmd with arguments len num arg
// len - length of the interval, (num + 1) - number of grid points
// ans - values of the function in the grid points
void calc_fun (string cmd, string arg, double len, int num, vector<double>& ans)
{
    string arg1 = toString(len) + " " + toString(num);
    string full_cmd = cmd + " " + arg1 + " " + arg;
    cout << full_cmd << endl;
    FILE* pf = popen(full_cmd.c_str(), "r");
    ans.resize(num + 1);
    for (int i = 0; i <= num; ++i)
        fscanf(pf, "%lf", &ans[i]);
    pclose(pf);
}

InputData::InputData (string input_file_name)
{
    try {
        po::options_description desc("Input data");
        desc.add_options()
            ("L", po::value<double>(), "Length of the space interval")
            ("T", po::value<double>(), "Length of the time interval")
            ("a", po::value<double>(), "Diffusion coefficient a")
            ("alpha", po::value<double>(), "Diffusion coefficient alpha")
            ("kappa_a", po::value<double>(), "Coefficient kappa_a")
            ("b", po::value<double>(), "Coefficient b")
            ("beta", po::value<double>(), "Coefficient beta (in BC)")
            ("gamma", po::value<double>(), "Coefficient gamma (in BC)")
            ("theta_b1", po::value<double>(), "Boundary temperature theta_b1")
            ("theta_b2", po::value<double>(), "Boundary temperature theta_b2")
            ("mode", po::value<string>(),
                "Mode: whether q(t) or r(t) is given (given_q or given_r)")
            ("q_init_guess", po::value<double>(),
                "Initial guess in the bisection method")
            ("q_init_len", po::value<double>(),
                "Initial length of the interval in the bisection method")
            ("q_tol", po::value<double>(), "Tolerance in the bisection method")
            ("verify_monotonicity", po::value<string>(),
                "Verify monotonicity (on or off)")
            ("monot_ver_q_1", po::value<double>(),
                "Left bound of the range of the monotonicity verification")
            ("monot_ver_q_2", po::value<double>(),
                "Right bound of the range of the monotonicity verification")
            ("monot_ver_q_step", po::value<double>(),
                "Step of the monotonicity verification")
            ("N", po::value<int>(), "Number of space grid points")
            ("M", po::value<int>(), "Number of time grid points")
            ("theta_0_fun_cmd", po::value<string>(), "Command theta_0")
            ("theta_0_fun_arg", po::value<string>(), "Arguments of theta_0")
            ("f_fun_cmd", po::value<string>(), "Command f")
            ("f_fun_arg", po::value<string>(), "Arguments of f")
            ("g_fun_cmd", po::value<string>(), "Command g")
            ("g_fun_arg", po::value<string>(), "Arguments of g")
            ("q_fun_cmd", po::value<string>(), "Command q")
            ("q_fun_arg", po::value<string>(), "Arguments of q")
            ("r_fun_cmd", po::value<string>(), "Command r")
            ("r_fun_arg", po::value<string>(), "Arguments of r")
            ("linear_sys_sol_method", po::value<string>(),
                "Linear system solution method (m - MTL or u - UMFPACK)")
            ("Newton_tol", po::value<double>(), "Tolerance in Newton's method")
            ("linear_sys_tol", po::value<double>(),
                "Tolerance in the linear system solution method")
        ;
        po::variables_map vm;
        ifstream ifs(input_file_name.c_str());
        if (!ifs)
            print_error("Can not open input file: " + input_file_name);
        else {
            po::store(parse_config_file(ifs, desc), vm);
            po::notify(vm);
        }
        get_pos_double_param(vm, "L", L);
        get_pos_double_param(vm, "T", T);
        get_pos_double_param(vm, "a", a);
        get_pos_double_param(vm, "alpha", alpha);
        get_pos_double_param(vm, "kappa_a", kappa_a);
        get_pos_double_param(vm, "b", b);
        get_pos_double_param(vm, "beta", beta);
        get_pos_double_param(vm, "gamma", gamma);
        get_pos_double_param(vm, "theta_b1", theta_b1);
        get_pos_double_param(vm, "theta_b2", theta_b2);
        string s;
        get_string_param(vm, "mode", s);
        if (s == "given_q")
            mode = MODE_GIVEN_Q;
        else if (s == "given_r")
            mode = MODE_GIVEN_R;
        else
            print_error("mode doesn't equal either given_q or given_r");
        get_double_param(vm, "q_init_guess", q_init_guess);
        get_pos_double_param(vm, "q_init_len", q_init_len);
        get_pos_double_param(vm, "q_tol", q_tol);
        get_string_param(vm, "verify_monotonicity", s);
        if (s == "on")
            verify_monotonicity = true;
        else if (s == "off")
            verify_monotonicity = false;
        else {
            print_error(
                "verify_monotonicity doesn't equal either on or off");
        }
        if (verify_monotonicity) {
            get_double_param(vm, "monot_ver_q_1", monot_ver_q_1);
            get_double_param(vm, "monot_ver_q_2", monot_ver_q_2);
            get_pos_double_param(vm, "monot_ver_q_step", monot_ver_q_step);
            if (monot_ver_q_1 >= monot_ver_q_2)
                print_error("monot_ver_q_1 >= monot_ver_q_2");
        }
        get_int_param(vm, "N", N);
        get_int_param(vm, "M", M);
        string theta_0_fun_cmd, theta_0_fun_arg, f_fun_cmd, f_fun_arg,
            g_fun_cmd, g_fun_arg, q_fun_cmd, q_fun_arg, r_fun_cmd, r_fun_arg;
        get_string_param(vm, "theta_0_fun_cmd", theta_0_fun_cmd);
        get_string_param(vm, "theta_0_fun_arg", theta_0_fun_arg);
        calc_fun(theta_0_fun_cmd, theta_0_fun_arg, L, N, theta_0);
        get_string_param(vm, "f_fun_cmd", f_fun_cmd);
        get_string_param(vm, "f_fun_arg", f_fun_arg);
        calc_fun(f_fun_cmd, f_fun_arg, L, N, f);
        get_string_param(vm, "g_fun_cmd", g_fun_cmd);
        get_string_param(vm, "g_fun_arg", g_fun_arg);
        calc_fun(g_fun_cmd, g_fun_arg, L, N, g);
        if (mode == MODE_GIVEN_Q) {
            get_string_param(vm, "q_fun_cmd", q_fun_cmd);
            get_string_param(vm, "q_fun_arg", q_fun_arg);
            calc_fun(q_fun_cmd, q_fun_arg, T, M, q);
        }
        else {  // mode == MODE_GIVEN_R
            get_string_param(vm, "r_fun_cmd", r_fun_cmd);
            get_string_param(vm, "r_fun_arg", r_fun_arg);
            calc_fun(r_fun_cmd, r_fun_arg, T, M, r);
        }
        get_pos_double_param(vm, "Newton_tol", Newton_tol);
        get_string_param(vm, "linear_sys_sol_method", s);
        if (s == "m")
            linear_sys_sol_method = SOL_METHOD_MTL;
        else if (s == "u")
            linear_sys_sol_method = SOL_METHOD_UMFPACK;
        else {
            print_error(
                "linear_sys_sol_method doesn't equal either m or u");
        }
        if (linear_sys_sol_method == SOL_METHOD_MTL)
            get_pos_double_param(vm, "linear_sys_tol", linear_sys_tol);
    }
    catch (ParamIsNotSet e) {
        print_error("Parameter \"" + e.param + "\" is not specified");
    }
    catch (ParamIsNotPositive e) {
        print_error("Parameter \"" + e.param + "\" is not positive");
    }
    catch (std::exception& e) {
        print_error(std::string("error: ") + e.what());
    }
    catch (...) {
        print_error("Exception of unknown type!");
    }
    // increase the upper bound a little bit
    // so that the for loop works correctly
    if (monot_ver_q_2 > 0)
        monot_ver_q_2 *= 1 + 1e-5;
    else
        monot_ver_q_2 *= 1 - 1e-5;
}
