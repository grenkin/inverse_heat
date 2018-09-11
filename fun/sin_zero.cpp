// Usage: sin_zero len num x1 x2 A
// A sin(2 pi (x - x1)/(x2 - x1)) if x in [x1, x2] and 0 else

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
    double x1 = atof(argv[3]);
    double x2 = atof(argv[4]);
    double A = atof(argv[5]);

    double h = len / num;
    for (int i = 0; i <= num; ++i) {
        double x = h * i;
        double y;
        if (x >= x1 && x <= x2)
            y = A * sin(2 * M_PI * (x - x1) / (x2 - x1));
        else
            y = 0.0;
        cout << y << " ";
    }

    return 0;
}
