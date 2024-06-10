#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <fstream>
#include <string>
#include <thread>
#include <sstream>
#include <algorithm>
#include <chrono>
#include <mutex>
#include <set>
#include <random>
#include <atomic>
#include "fft.cpp"

typedef std::complex<double> complex;

template<typename T>
T evaluate_polynomial_trivial(std::vector<T> P, T x){
    int N = P.size();
    T X = 1;
    T y;
    for (int j = 0; j < N; j++){
        y += P[j] * X;
        X *= x;
    }
    return y;
}

template<typename T>
std::vector<T> add_polynomials(std::vector<T> P, std::vector<T> Q){
    int N = std::max(P.size(), Q.size());
    P.resize(N);
    Q.resize(N);
    std::vector<T> R(N);
    for (int j = 0; j < N; j++){
        R[j] = P[j] + Q[j];
    }
    return R;
}

template<typename T>
std::vector<T> multiply_polynomials_trivial(std::vector<T> P, std::vector<T> Q){
    int N = P.size();
    int M = Q.size();
    std::vector<T> R(M + N, 0.);
    for (int j = 0; j < N; j++){
        for (int k = 0; k < M; k++){
            R[j + k] += P[j] * Q[k];
        }
    }
    return R;
}

template<typename T>
T horner_evaluate(std::vector<T> P, T x){
    int N = P.size();
    T y = P[N - 1];
    for (int j = N - 2; j >= 0; j--){
        y = x * y + P[j];
    }
    return y;
}



std::vector<complex> multiply_polynomials_FFT(std::vector<complex> P, std::vector<complex> Q){
    int N = P.size();
    int M = Q.size();
    int l = next_power_of_two(M + N);
    P.resize(l);
    Q.resize(l);
    std::vector<complex> P_star = Radix2FFT(P);
    std::vector<complex> Q_star = Radix2FFT(Q);
    std::vector<complex> R_star(l);
    for (int j = 0; j < l; j++){
        R_star[j] = P_star[j] * Q_star[j];
    }
    return InverseRadix2FFT(R_star);
}