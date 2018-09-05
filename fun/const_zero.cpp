// Usage: const_zero len num x1 x2 val

#include <iostream>
#include <cstdlib>

using namespace std;

int main (int argc, char **argv)
{
    double len = atof(argv[1]);
    int num = atoi(argv[2]);
    double x1 = atof(argv[3]);
    double x2 = atof(argv[4]);
    double val = atof(argv[5]);

    double h = len / num;
    for (int i = 0; i <= num; ++i) {
        double x = h * i;
        if (x >= x1 && x <= x2)
            cout << val;
        else
            cout << 0.0;
        cout << " ";
    }

    return 0;
}
