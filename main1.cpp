#include <iostream> // for console output
#include <random> // for random number generator
#include <fstream> // for files output
using namespace std;
#define e 7.5 // edge of the box: (0; e)x(0; e)x(0; e)
#define n 125 // number of particles
#define n3 5
#define dt 0.001 // dt
#define v0 2.55 // initial velocity
#define N 1000 // number of iterations
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
    double er = 1 / e;
    double nr = 1 / n;
    double msd = 0;
    ofstream of("pos.xyz"); // creating a file
    ofstream af("vel.csv");
    ofstream ef("eng.csv");
    ofstream mf("msd.csv");
    mf << "N,msd" << endl << "0,0" << endl;
    ef << "F,K,P" << endl;
    af << "vx,vy,vz" << endl;
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
    of << n << endl << endl;
    for(int i = 0; i < n3; i++) {
        for(int j = 0; j < n3; j++) {
            for(int k = 0; k < n3; k++) {
                int t = k + j * n3 + i * n3 * n3;
                vec vtmp(0.5+i, 0.5+j, 0.5+k);
                R0[t] = vtmp*e/n3;
                R00[t] = R0[t];
                R000[t] = R0[t];
                of << "A" << t << " " << R0[t].x << " " << R0[t].y << " " << R0[t].z << endl;
            }
        }
    }
    for(int i = 0; i < n; i++) {
        double a = dis(gen);
        double b = dis(gen);
        vec vtmp(cos(b) * cos(a), cos(b) * sin(a), sin(b));
        V[i] = v0 * vtmp;
        vtemp = vtemp + V[i]*nr;
    }
    for(int i = 0; i < n; i++) {
        V[i] = V[i] - vtemp;
        for(int j = 0; j < n; j++) {
            if(i != j) {
                vec v = R0[i] - R0[j];
                vec diff(round(v.x * er), round(v.y * er), round(v.z * er));
                v = v - e * diff;
                double r2 = mod2(v);
                A0[i] = A0[i] + 24 * v * (2 - pow(r2, 3)) / pow(r2, 7);  
            }
        }
    }
    // calculating R1
    of << n << endl << endl;
    for(int i = 0; i < n; i++) {
        double a = dis(gen);
        double b = dis(gen);
        R1[i] = R0[i] + V[i] * dt + 0.5 * A0[i] * dt * dt;
        R11[i] = R1[i];
        msd += mod2(R1[i]-R000[i])*nr;
        of << "A" << i << " " << R1[i].x << " " << R1[i].y << " " << R1[i].z << endl;
    }
    mf << 1 << "," << msd << endl;
    //
    for(int I = 0; I < N; I++) { // iterating
        of << n << endl << I << endl;
        for(int i = 0; i < n; i++) {
            vec zero(0, 0, 0);
            A[i] = zero;
            for(int j = 0; j < n; j++) {
                if(i != j) {
                    vec v = R1[i] - R1[j];
                    vec diff(round(v.x * er), round(v.y * er), round(v.z * er));
                    v = v - e * diff;
                    double r2 = mod2(v);
                    A[i] = A[i] + 24 * v * (2 - pow(r2, 3)) / pow(r2, 7);
                }
            }
        }
        msd = 0;
        for(int i = 0; i < n; i++) { // updating R1, V, R0
            of << "A" << i << " " << R1[i].x << " " << R1[i].y << " " << R1[i].z << endl;
            vec Rt = R1[i];
            vec Rtt = R11[i];
            R1[i] = pos(R0[i], R1[i], A[i]);
            R11[i] = pos(R00[i], R11[i], A[i]);
            msd += mod2(R11[i]-R000[i])*nr;
            V[i] = vel(R0[i], R1[i]);
            R0[i] = Rt;
            R00[i] = Rtt;
            // if particle escapes then move it to the other side of the box
            vec diff(floor(R1[i].x * er), floor(R1[i].y * er), floor(R1[i].z * er));
            R0[i] = R0[i] - e * diff;
            R1[i] = R1[i] - e * diff;
        }
        mf << I+2 << "," << msd << endl;
        K = 0;
        P = 0;
        for(int i = 0; i < n; i++) { // calculating energy
            K += mod2(V[i]);
            for(int j = i+1; j < n; j++) {
                if(i != j) {
                    vec v = R0[i] - R0[j];
                    vec diff(round(v.x * er), round(v.y * er), round(v.z * er));
                    v = v - e * diff;
                    double r2 = mod2(v);
                    P += 4 * (1 - pow(r2, 3)) / pow(r2, 6);
                }
            }
        }
        K *= 0.5;
        ef << K+P << "," << K << "," << P << endl;
    }
    for(int i = 0; i < n; i++) {
        af << V[i] << endl;
    }
    of.close(); // closing the file
    af.close();
    ef.close();
    return 0;
}