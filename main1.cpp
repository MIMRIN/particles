#include <iostream> // for console output
#include <random> // for random number generator
#include <fstream> // for files output
using namespace std;
#define e 9 // edge of the box: (0; e)x(0; e)x(0; e)
#define c 3 // cut-off radius
#define n 216 // number of particles
#define n3 6
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
        os << num.x << " " << num.y << " " << num.z;
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
    ofstream of("list.xyz"); // creating a file
    double K = 0; // kinetic energy
    double P = 0; // potential energy
    vec R0[n]; // initial positions
    vec A0[n]; // initial accelerations
    vec R1[n]; // 1-step positions
    vec V[n]; // 1 step velocities
    // random number generator
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dis(0, 2 * M_PI); // from 0 to 2pi
    // generating initial positions
    of << n << endl << endl;
    for(int i = 1; i <= n3; i++) {
            for(int j = 1; j <= n3; j++) {
                for(int k = 1; k <= n3; k++) {
                    int t = k + (j-1) * n3 + (i-1) * n3 * n3 - 1;
                    R0[t].x = e / (n3+1) * i;
                    R0[t].y = e / (n3+1) * j;
                    R0[t].z = e / (n3+1) * k;
                    of << "A" << t << " " << R0[t].x << " " << R0[t].y << " " << R0[t].z << endl;
                }
            }
    }    
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            if(i != j) {
                vec v = R0[i] - R0[j];
                v.x -= e * round(v.x / e);
                v.y -= e * round(v.y / e);
                v.z -= e * round(v.z / e);
                double r2 = mod2(v);
                if(r2 <= c * c) {
                    A0[i] = A0[i] + 24 * v * (2 - pow(r2, 3)) / pow(r2, 7);
                }
            }
        }
    }
    // calculating R1
    of << n << endl << endl;
    for(int i = 0; i < n; i++) {
        double a = dis(gen);
        double b = dis(gen);
        R1[i].x = R0[i].x + v0 * cos(b) * cos(a) * dt + 0.5 * dt * dt * A0[i].x;
        R1[i].y = R0[i].y + v0 * cos(b) * sin(a) * dt + 0.5 * dt * dt * A0[i].y;
        R1[i].z = R0[i].z + v0 * sin(b) * dt + 0.5 * dt * dt * A0[i].z;
        of << "A" << i << " " << R1[i].x << " " << R1[i].y << " " << R1[i].z << endl;
    }
    //
    for(int I = 0; I < N; I++) { // iterating
        of << n << endl << I << endl;
        vec A[n]; // initiating acceleration
        for(int i = 0; i < n; i++) {
            for(int j = 0; j < n; j++) {
                if(i != j) {
                    vec v = R1[i] - R1[j];
                    v.x -= e * round(v.x / e);
                    v.y -= e * round(v.y / e);
                    v.z -= e * round(v.z / e);
                    double r2 = mod2(v);
                    if(r2 <= c * c) {
                        A[i] = A[i] + 24 * v * (2 - pow(r2, 3)) / pow(r2, 7);
                    }
                }
            }
        }
        for(int i = 0; i < n; i++) { // updating R1, V, R0
            of << "A" << i << " " << R1[i].x << " " << R1[i].y << " " << R1[i].z << endl;
            vec Rt = R1[i];
            R1[i] = pos(R0[i], R1[i], A[i]);
            V[i] = vel(R0[i], R1[i]);
            R0[i] = Rt;
            // if particle escapes then move it to the other side of the box
            R0[i].x -= e * floor(R1[i].x / e);
            R1[i].x -= e * floor(R1[i].x / e);
            R0[i].y -= e * floor(R1[i].y / e);
            R1[i].y -= e * floor(R1[i].y / e);
            R0[i].z -= e * floor(R1[i].z / e);
            R1[i].z -= e * floor(R1[i].z / e);

        }         
        K = 0;
        P = 0;
        for(int i = 0; i < n; i++) { // calculating energy
            K += 0.5 * mod2(V[i]);
            for(int j = i+1; j < n; j++) {
                vec v = R0[i] - R0[j];
                v.x -= e * round(v.x / e);
                v.y -= e * round(v.y / e);
                v.z -= e * round(v.z / e);
                double r2 = mod2(v);
                P += 4 * (1 - pow(r2, 3)) / pow(r2, 6);
            }
        }
        cout << I << " " << K + P << endl;
    }
    of.close(); // closing the file
    return 0;
}