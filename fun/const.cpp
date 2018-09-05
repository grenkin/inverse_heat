// Usage: const len num val

#include <iostream>
#include <cstdlib>

using namespace std;

int main (int argc, char **argv)
{
    double len = atof(argv[1]);
    int num = atoi(argv[2]);
    double val = atof(argv[3]);

    for (int i = 0; i <= num; ++i)
        cout << val << " ";

    return 0;
}
