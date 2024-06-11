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

void compute_multiplication_segment(const std::vector<complex>& P, const std::vector<complex>& Q, std::vector<complex>& R, size_t start_index, size_t num_threads, size_t n){
    for (size_t k = start_index; k < n ; k += num_threads) {
        R[k] = P[k] * R[k];
    }
}


std::vector<complex> multiply_polynomials_FFT_parallel(std::vector<complex> P, std::vector<complex> Q, int num_threads){
    int N = P.size();
    int M = Q.size();
    int l = next_power_of_two(M + N);
    P.resize(l);
    Q.resize(l);
    FFT_Radix2_Parallel(P, num_threads / 2);
    FFT_Radix2_Parallel(Q, num_threads - num_threads / 2);
    std::vector<complex> R(l);

    std::vector<std::thread> threads;
    threads.reserve(num_threads - 1);

    for (size_t i = 1; i < num_threads; ++i) {
        threads.emplace_back(compute_multiplication_segment, std::cref(P), std::cref(Q), std::ref(R), i, num_threads, l);
    }
    compute_multiplication_segment(P, Q, R, 0, num_threads, l);

    for (auto& th : threads) {
        th.join();
    }
    for (int j = 0; j < l; j++){
        R[j] = P[j] * Q[j];
    }
    Inverse_Radix2_Parallel(R, num_threads);
    return R;
}