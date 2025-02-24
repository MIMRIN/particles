#include <iostream>
#include <cmath>
#include <random>
#include <fstream>
using namespace std;
#define e 7 // edge of the box
#define c 3 // cut-off radius
#define n 100 // number of particles
#define dt 0.01 // dt
#define v0 2 // initial velocity
#define N 100
// defining vectors
class vec {
public:
    long double x, y, z;
    // Конструктор
    vec() : x(0), y(0), z(0) {}
    vec(long double X, long double Y, long double Z) : x(X), y(Y), z(Z) {}
    // Перегрузка оператора сложения
    vec operator+(const vec& other) const {
        return vec(x + other.x, y + other.y, z + other.z);
    }

    // Перегрузка оператора вычитания
    vec operator-(const vec& other) const {
        return vec(x - other.x, y - other.y, z - other.z);
    }

    // Перегрузка оператора умножения
    vec operator*(const long double& other) const {
        return vec(x * other, y * other, z * other);
    }
    // Перегрузка оператора деления
    vec operator/(const long double& other) const {
        return vec(x / other, y / other, z / other);
    }
    // Перегрузка оператора вывода для удобства
    friend std::ostream& operator<<(std::ostream& os, const vec& num) {
        os << num.x << " " << num.y << " " << num.z;
        return os;
    }
};
vec operator*(long double scalar, const vec& v) {
    return v * scalar; // Используем уже перегруженный оператор vec * long double
}
// defining vector functions
vec pos(vec r0, vec r1, vec a1) {
    return 2*r1-r0+a1*dt*dt;
}
long double mod(vec v) {
    return pow(v.x*v.x + v.y*v.y + v.z*v.z, 0.5);
}
vec vel(vec r1, vec r3) {
    return (r3 - r1) / (2 * dt);
}

int main()
{
    ofstream of("list.xyz");
    long double K = 0;
    long double P = 0;
    // creating an array for PBC
    vec rrr[26];
    int kkk = 0;
    for(int i1 = -1; i1 <= 1; i1++) {
        for(int i2 = -1; i2 <= 1; i2++) {
            for(int i3 = -1; i3 <= 1; i3++) {
                if(i1 != 0 || i2 != 0 || i3 != 0) {
                    vec rr(i1 * e, i2 * e, i3 * e);
                    rrr[kkk] = rr;
                    kkk++;
                }
            }
        }
    }
    //
    vec R0[n]; // initial positions
    vec A0[n]; // initial accelerations
    vec R1[n]; // 1-step positions
    vec V[n]; // 1 step velocities
    //
    // random number generator
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(0, e);
    // generating initial positions
    of << n << endl << endl;
    for(int i = 0; i < n; i++) {
        R0[i].x = dis(gen);
        R0[i].y = dis(gen);
        R0[i].z = dis(gen);
        of << "A" << i << " " << R0[i].x << " " << R0[i].y << " " << R0[i].z << endl;
    }
    //
    for(int i = 0; i < n; i++) { // calculating A0
        for(int j = 0; j < n && j != i; j++) {
            vec v = R0[i] - R0[j];
            long double r = mod(v);
            if(r <= c) {
                A0[i] = A0[i] + v * (2 - pow(r, 6)) / pow(r, 14);
            }
            else {
                for(int k = 0; k < 26; k++) {
                    vec v1 = rrr[k] + v;
                    long double r1 = mod(v1);
                    if(r1 <= c) {
                        A0[i] = A0[i] + v1 * (2 - pow(r1, 6)) / pow(r1, 14);
                        break; // if found then no need to continue
                    }
                }
            }
        }
        A0[i] = A0[i] * 24;
    }
    // находим R1
    of << n << endl << endl;
    for(int i = 0; i < n; i++) {
        long double a = dis(gen) / e * 2 * M_PI;
        long double b = dis(gen) / e * 2 * M_PI;
        R1[i].x = R0[i].x + v0 * cos(a) * cos(b) * dt + 0.5 * dt * dt * A0[i].x;
        R1[i].y = R0[i].y + v0 * cos(b) * sin(a) * dt + 0.5 * dt * dt * A0[i].y;
        R1[i].z = R0[i].z + v0 * sin(b) * dt + 0.5 * dt * dt * A0[i].z;
        of << "A" << i << " " << R1[i].x << " " << R1[i].y << " " << R1[i].z << endl;
    }
    //
    for(int I = 0; I < N; I++) { // совершаем итерации
        of << n << endl << endl;
        vec A[n]; // инициализируем нулевое ускорение
        for(int i = 0; i < n; i++) { // находим ускорение для следрадвека
            for(int j = 0; j < n && j != i; j++) {
                vec v = R1[i] - R1[j];
                long double r = mod(v);
                if(r <= c) {
                    A[i] = A[i] + v * (2 - pow(r, 6)) / pow(r, 14);
                }
                else {
                    for(int k = 0; k < 26; k++) {
                        vec v1 = rrr[k] + v;
                        long double r1 = mod(v1);
                        if(r1 <= c) {
                            A[i] = A[i] + v1 * (2 - pow(r1, 6)) / pow(r1, 14);
                            break; // если нашли, то дальше не ищем
                        }
                    }
                }
            }
            A[i] = A[i] * 24;
        }
        for(int i = 0; i < n; i++) { // обновляем радвеки
            of << "A" << i << " " << R1[i].x << " " << R1[i].y << " " << R1[i].z << endl;
            vec Rt = R1[i];
            R1[i] = pos(R0[i], R1[i], A[i]);
            V[i] = vel(R0[i], R1[i]);
            R0[i] = Rt;
            // если вышел из коробки, то корректируем данные
            if(R1[i].x > e) {
                R0[i].x = R0[i].x - R1[i].x;
                R1[i].x = 0;
            }
            else if(R1[i].x < 0) {
                R0[i].x = R0[i].x - R1[i].x + e;
                R1[i].x = e;
            }
            if(R1[i].y > e) {
                R0[i].y = R0[i].y - R1[i].y;
                R1[i].y = 0;
            }
            else if(R1[i].y < 0) {
                R0[i].y = R0[i].y - R1[i].y + e;
                R1[i].y = e;
            }    
            if(R1[i].z > e) {
                R0[i].z = R0[i].z - R1[i].z;
                R1[i].z = 0;
            }
            else if(R1[i].z < 0) {
                R0[i].z = R0[i].z - R1[i].z + e;
                R1[i].z = e;
            }
        }
        K = 0;
        P = 0;
        for(int i = 0; i < n; i++) { //считаем полную энергию
            K += pow(mod(V[i]), 2);
            for(int j = i+1; j < n; j++) {
                long double r = mod(R0[i] - R0[j]);
                P += (1 - pow(r, 6)) / pow(r, 12);
            }
        }
        K *= 0.5;
        P *= 4;
        cout << K << endl;
    }
    of.close();
    return 0;
}