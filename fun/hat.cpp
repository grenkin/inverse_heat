// Usage: hat len num x1 x2 val

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
    double val = atof(argv[5]);

    double x_mid = (x1 + x2) / 2;
    double h = len / num;
    for (int i = 0; i <= num; ++i) {
        double x = h * i;
        double y;
        if (x < x1 || x > x2)
            y = 0.0;
        else if (x < x_mid)
            y = val * (x - x1) / (x_mid - x1);
        else
            y = val * (x - x2) / (x_mid - x2);
        cout << y << " ";
    }

    return 0;
}
