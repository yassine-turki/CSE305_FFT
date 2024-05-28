#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#define M_PI 3.14159265358979323846

typedef std::complex<double> complex;

std::vector<complex> Radix2FFT(std::vector<complex> P) {
    //Cooley-Tukey Radix-2 FFT

    int N = P.size();
    if (N == 1) {
        return P;
    }
    std::vector<complex> U(N / 2), V(N / 2);
    for (int i = 0; i < N / 2; i++) {
        U[i] = P[2 * i];
        V[i] = P[2 * i + 1];
    }

    std::vector<complex> U_star = Radix2FFT(U);
    std::vector<complex> V_star = Radix2FFT(V);
    std::vector<complex> res(N);
    double theta = (2 * M_PI) / N;
    complex omega_n = std::polar(1.0, theta);
    complex omega = 1;
    for (int i = 0; i < N / 2; i++) {
        res[i] = U_star[i] + omega * V_star[i];
        res[i + N / 2] = U_star[i] - omega * V_star[i];
        omega *= omega_n;
    }
    return res;
}

std::vector<complex> InverseRadix2FFT(std::vector<complex> P_star) {
    int N = P_star.size();
    if (N == 1){
        return P_star;
    }

    std::vector<complex> U_star;
    std::vector<complex> V_star;

    for (int j = 0; j < N / 2; j++){
        U_star[j] = P_star[2 * j];
        V_star[j] = P_star[2 * j + 1];
    }

    std::vector<complex> U = InverseRadix2FFT(U_star);
    std::vector<complex> V = InverseRadix2FFT(V_star);

    double theta = 2 * M_PI - (2 * M_PI) / N;
    complex omega_n = std::polar(1.0, theta);
    complex omega = 1;
    std::vector<complex> P(N); 

    for (int j = 0; j < N / 2; j++){
        P[j] = (U[j] + omega * V[j]) / 2.;
        P[j + N / 2] = (U[j] - omega * V[j]) / 2.;
        omega = omega * omega_n;
    }

    return P;
}

int main() {
    std::vector<complex> input = {1.0, 2.0, 3.0, 4.0};
    std::vector<complex> result = Radix2FFT(input);
    for (const auto& val : result) {
        std::cout << val << " ";
    }
    std::cout << std::endl;
    return 0;
}