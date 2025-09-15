#include <iostream> // for console output
#include <random> // for random number generator
#include <fstream> // for files output
using namespace std;
#define e 4.5 // edge of the box: (0; e)x(0; e)x(0; e)
#define n 27 // number of particles
#define n3 3
#define dt 0.001 // dt
#define v0 2 // initial velocity
#define N 10000 // number of iterations
// defining vectors
class vec {
public:
    double x, y, z;
    // Конструктор
    vec() : x(0), y(0), z(0) {}
    vec(double X, double Y, double Z) : x(X), y(Y), z(Z) {}
    // Перегрузка оператора сложения
    vec operator+(const vec& other) const {
        return vec(x + other.x, y + other.y, z + other.z);
    }

    // Перегрузка оператора вычитания
    vec operator-(const vec& other) const {
        return vec(x - other.x, y - other.y, z - other.z);
    }

    // Перегрузка оператора умножения
    vec operator*(const double& other) const {
        return vec(x * other, y * other, z * other);
    }
    // Перегрузка оператора деления
    vec operator/(const double& other) const {
        return vec(x / other, y / other, z / other);
    }
    // Перегрузка оператора вывода для удобства
    friend std::ostream& operator<<(std::ostream& os, const vec& num) {
        os << num.x << "," << num.y << "," << num.z;
        return os;
    }
};
vec operator*(double scalar, const vec& v) {
    return v * scalar; // Используем уже перегруженный оператор vec * double
}
// defining vector functions
vec pos(vec r0, vec r1, vec a1) { // verlet integration
    return 2*r1-r0+a1*dt*dt;
}
double mod2(vec v) { // module
    return v.x*v.x + v.y*v.y + v.z*v.z;
}
vec vel(vec r0, vec r2) { // verlet velocities
    return (r2 - r0) / (2 * dt);
}

int main() {
    ofstream f1("animation.xyz");
    ofstream f2("data.csv");
    f2 << "x,y,z,X,Y,Z\n";
    double K = 0; // kinetic energy
    double P = 0; // potential energy
    vec vtemp(0, 0, 0);
    vec R0[n]; // initial positions
    vec R00[n];
    vec R000[n];
    vec R11[n];
    vec A0[n]; // initial accelerations
    vec A[n];
    vec R1[n]; // 1-step positions
    vec V[n]; // 1 step velocities
    // random number generator
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(0, 2 * M_PI); // from 0 to 2pi
    // generating initial positions
    f1 << n << "\n\n";
    for(int i = 0; i < n3; i++) {
        for(int j = 0; j < n3; j++) {
            for(int k = 0; k < n3; k++) {
                int t = k + j * n3 + i * n3 * n3;
                vec vtmp(0.5+i, 0.5+j, 0.5+k);
                R0[t] = vtmp*e/n3;
                R00[t] = R0[t];
                R000[t] = R0[t];
                f1 << "A" << t << " " << R0[t].x << " " << R0[t].y << " " << R0[t].z << "\n";
                f2 << R0[t] << "," << R00[t] << "\n";
            }
        }
    }
    for(int i = 0; i < n; i++) {
        double a = dis(gen);
        double b = dis(gen);
        vec vtmp(cos(b) * cos(a), cos(b) * sin(a), sin(b));
        V[i] = v0 * vtmp;
        vtemp = vtemp + V[i]/n;
    }
    for(int i = 0; i < n; i++) {
        V[i] = V[i] - vtemp;
        for(int j = 0; j < n; j++) {
            if(i != j) {
                vec v = R0[i] - R0[j];
                vec diff(round(v.x/e), round(v.y/e), round(v.z/e));
                v = v - e * diff;
                double r2 = mod2(v);
                A0[i] = A0[i] + 24 * v * (2 - pow(r2, 3)) / pow(r2, 7);  
            }
        }
    }
    // calculating R1
    f1 << n << "\n\n";
    for(int i = 0; i < n; i++) {
        R1[i] = R0[i] + V[i] * dt + 0.5 * A0[i] * dt * dt;
        R11[i] = R1[i];
        f1 << "A" << i << " " << R1[i].x << " " << R1[i].y << " " << R1[i].z << "\n";
        f2 << R1[i] << "," << R11[i] << "\n";
    }
    //
    for(int I = 2; I < N; I++) { // iterating
        f1 << n << "\n\n";
        for(int i = 0; i < n; i++) {
            vec zero(0, 0, 0);
            A[i] = zero;
            for(int j = 0; j < n; j++) {
                if(i != j) {
                    vec v = R1[i] - R1[j];
                    vec diff(round(v.x/e), round(v.y/e), round(v.z/e));
                    v = v - e * diff;
                    double r2 = mod2(v);
                    A[i] = A[i] + 24 * v * (2 - pow(r2, 3)) / pow(r2, 7);
                }
            }
        }
        for(int i = 0; i < n; i++) { // updating R1, V, R0
            vec Rt = R1[i];
            vec Rtt = R11[i];
            R1[i] = pos(R0[i], R1[i], A[i]);
            R11[i] = pos(R00[i], R11[i], A[i]);
            R0[i] = Rt;
            R00[i] = Rtt;
            // if particle escapes then move it to the other side of the box
            vec diff(floor(R1[i].x/e), floor(R1[i].y/e), floor(R1[i].z/e));
            R0[i] = R0[i] - e * diff;
            R1[i] = R1[i] - e * diff;
            f1 << "A" << i << " " << R1[i].x << " " << R1[i].y << " " << R1[i].z << "\n";
            f2 << R1[i] << "," << R11[i] << "\n";
        }
    }
    return 0;
}