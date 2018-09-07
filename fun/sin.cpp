// Usage: sin len num A k phi
// A sin(k pi/L x + phi)

#include <cmath>
#include <iostream>
#include <cstdlib>

using namespace std;

int main (int argc, char **argv)
{
    if (argc != 6) {
        cerr << "Wrong arguments!" << endl;
        exit(1);
    }
    double len = atof(argv[1]);
    int num = atoi(argv[2]);
    double A = atof(argv[3]);
    double k = atof(argv[4]);
    double phi = atof(argv[5]);

    double h = len / num;
    for (int i = 0; i <= num; ++i) {
        double x = h * i;
        double y = A * sin(k * M_PI / len * x + phi);
        cout << y << " ";
    }

    return 0;
}
