// Usage: const len num val

#include <iostream>
#include <cstdlib>

using namespace std;

int main (int argc, char **argv)
{
    if (argc != 4) {
        cerr << "Wrong arguments!" << endl;
        exit(1);
    }
    double len = atof(argv[1]);
    int num = atoi(argv[2]);
    double val = atof(argv[3]);

    for (int i = 0; i <= num; ++i)
        cout << val << " ";

    return 0;
}
