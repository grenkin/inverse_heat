// Usage: sin_const len num C A k phi
// C + A sin(k pi/L x + phi)

#include <cmath>
#include <iostream>
#include <cstdlib>

using namespace std;

int main (int argc, char **argv)
{
    if (argc != 7) {
        cerr << "Wrong arguments!" << endl;
        exit(1);
    }
    double len = atof(argv[1]);
    int num = atoi(argv[2]);
    double C = atof(argv[3]);
    double A = atof(argv[4]);
    double k = atof(argv[5]);
    double phi = atof(argv[6]);

    double h = len / num;
    for (int i = 0; i <= num; ++i) {
        double x = h * i;
        double y = C + A * sin(k * M_PI / len * x + phi);
        cout << y << " ";
    }

    return 0;
}
