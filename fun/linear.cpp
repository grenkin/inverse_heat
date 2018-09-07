// Usage: linear len num val1 val2

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
    double val1 = atof(argv[3]);
    double val2 = atof(argv[4]);

    double h = len / num;
    for (int i = 0; i <= num; ++i) {
        double x = h * i;
        double val = val1 + x / len * (val2 - val1);
        cout << val << " ";
    }

    return 0;
}
