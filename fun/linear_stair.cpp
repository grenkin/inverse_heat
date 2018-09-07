// Usage: linear_stair len num x1 x2 val1 val2

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
    double x1 = atof(argv[3]);
    double x2 = atof(argv[4]);
    double val1 = atof(argv[5]);
    double val2 = atof(argv[6]);

    double h = len / num;
    for (int i = 0; i <= num; ++i) {
        double x = h * i;
        double y;
        if (x <= x1)
            y = val1;
        else if (x >= x2)
            y = val2;
        else
            y = val1 * (x2 - x) / (x2 - x1) + val2 * (x - x1) / (x2 - x1);
        cout << y << " ";
    }

    return 0;
}
