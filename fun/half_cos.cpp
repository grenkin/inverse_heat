// Usage: half_cos len num theta_b1 theta_b2
// (theta_b1 + theta_b2) / 2 - (theta_b2 - theta_b1) / 2 * cos(pi x / L)

#include <cmath>
#include <iostream>
#include <cstdlib>

using namespace std;

int main (int argc, char **argv)
{
    if (argc != 5) {
        cerr << "Wrong arguments!" << endl;
        exit(1);
    }
    double len = atof(argv[1]);
    int num = atoi(argv[2]);
    double theta_b1 = atof(argv[3]);
    double theta_b2 = atof(argv[4]);

    double h = len / num;
    for (int i = 0; i <= num; ++i) {
        double x = h * i;
        double y = (theta_b1 + theta_b2) / 2
            - (theta_b2 - theta_b1) / 2 * cos(M_PI * x / len);
        cout << y << " ";
    }

    return 0;
}
