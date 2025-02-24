#include <iostream>
#include <cstdint>
#include <climits>
using namespace std;

unsigned long long align(unsigned long long U, int P) {
    int p = 1;
    for(int i = 0; i < P; i++) {
        p *= 2;
    }
    if(p >= U) {
        return p;
    }
    else {
        if(U == ULLONG_MAX) {
            return ULLONG_MAX;
        }
        else if(U % p == 0) {
            return U;
        }
        else {
            return (U / p) * (p+1);
        }
    }
}
int main()
{
    cout << align(18446744073709551615, 16);
    return 0;
}