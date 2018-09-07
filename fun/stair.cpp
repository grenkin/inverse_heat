// Usage: stair len num x val1 val2

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
    double val1 = atof(argv[4]);
    double val2 = atof(argv[5]);

    double h = len / num;
    for (int i = 0; i <= num; ++i) {
        double x = h * i;
        double y;
        if (x <= x1)
            y = val1;
        else
            y = val2;
        cout << y << " ";
    }

    return 0;
}
