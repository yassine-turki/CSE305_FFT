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


#define M_PI 3.14159265358979323846

typedef std::complex<double> complex;



bool is_power_of_two(int n){
    return (n & (n - 1)) == 0;
}

int next_power_of_two(int n){
    int m = std::ceil(std::log2(n));
    return std::pow(2, m);
}


std::vector<complex> Naive_DFT(std::vector<complex> P) {
    int N = P.size();
    std::vector<complex> res(N);
    for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < N; j++) {
            complex t = std::polar(1.0, -2 * M_PI * i * j / N);
            res[i] += P[j] * t;
        }
    }
    return res;
}

std::vector<complex> Naive_IDFT(const std::vector<complex>& P) {
    int N = P.size();
    std::vector<complex> res(N);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            complex t = std::polar(1.0, 2 * M_PI * i * j / N);
            res[i] += P[j] * t;
        }
        res[i] /= N;
    }
    return res;
}

void compute_dft_segment(const std::vector<complex>& P, std::vector<complex>& res, size_t start_index, size_t num_threads, int N, std::mutex& mtx) {
    for (size_t k = start_index; k < N; k += num_threads) {
        complex sum = 0;
        for (size_t j = 0; j < N; ++j) {
            complex t = std::polar(1.0, -2 * M_PI * j * k / N);
            sum += P[j] * t;
        }
        std::lock_guard<std::mutex> lock(mtx);
        res[k] = sum;
    }
}


std::vector<complex> Naive_DFT_Parallel(std::vector<complex> &P, size_t num_threads) {
    int N = P.size();
    if (num_threads == 1) {
        return Naive_DFT(P);
    }

    std::vector<std::thread> threads;
    threads.reserve(num_threads - 1);

    std::vector<complex> res(N);
    std::mutex mtx;

    for (size_t i = 1; i < num_threads; ++i) {
        threads.emplace_back(compute_dft_segment, std::cref(P), std::ref(res), i, num_threads, N, std::ref(mtx));
    }
    compute_dft_segment(P, res, 0, num_threads, N, mtx);

    for (auto& th : threads) {
        th.join();
    }

    return res;
}

void compute_idft_segment(const std::vector<complex>& P, std::vector<complex>& res, size_t start_index, size_t num_threads, int N, std::mutex& mtx) {
    for (size_t k = start_index; k < N; k += num_threads) {
        complex sum = 0;
        for (size_t j = 0; j < N; ++j) {
            complex t = std::polar(1.0, 2 * M_PI * j * k / N);
            sum += P[j] * t;
        }
        std::lock_guard<std::mutex> lock(mtx);
        res[k] = sum / (double) N;
    }
}

std::vector<complex> Naive_IDFT_Parallel(std::vector<complex>& P, size_t num_threads) {
    int N = P.size();
    if (num_threads == 1) {
        return Naive_IDFT(P);
    }

    std::vector<std::thread> threads;
    threads.reserve(num_threads - 1);

    std::vector<complex> res(N);
    std::mutex mtx;

    for (size_t i = 1; i < num_threads; ++i) {
        threads.emplace_back(compute_idft_segment, std::cref(P), std::ref(res), i, num_threads, N, std::ref(mtx));
    }
    compute_idft_segment(P, res, 0, num_threads, N, mtx);

    for (auto& th : threads) {
        th.join();
    }

    return res;
}


std::vector<complex> Radix2FFT(std::vector<complex> P) {
    //Cooley-Tukey Radix-2 FFT
    int N = P.size();
    if(!is_power_of_two(N)){
        std::cout << "The input size " << N << " is not a power of 2! Resizing input" << std::endl;
        P.resize(next_power_of_two(P.size()));
        N = P.size();
    }

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
    double theta = (-2 * M_PI) / N;
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
    if (N == 1) {
        return P_star;
    }

    if (!is_power_of_two(N)) {
        std::cout << "The input size " << N << " is not a power of 2! Resizing input" << std::endl;
        P_star.resize(next_power_of_two(N));
        N = P_star.size();
    }

    std::vector<complex> U_star(N / 2), V_star(N / 2);

    for (int j = 0; j < N / 2; j++) {
        U_star[j] = P_star[2 * j];
        V_star[j] = P_star[2 * j + 1];
    }

    std::vector<complex> U = InverseRadix2FFT(U_star);
    std::vector<complex> V = InverseRadix2FFT(V_star);

    double theta = (-2 * M_PI) / N;
    complex omega_n = std::polar(1.0, theta);
    omega_n = std::conj(omega_n);
    complex omega = 1.0;
    std::vector<complex> P(N);

    for (int j = 0; j < N / 2; j++) {
        P[j] = (U[j] + omega * V[j]) / 2.0;
        P[j + N / 2] = (U[j] - omega * V[j]) / 2.0;
        omega *= omega_n;
    }

    return P;
}

void split_data_FFT(const std::vector<complex>& P, std::vector<complex>& U, std::vector<complex>& V, size_t start, size_t end) {
    for (size_t i = start; i < end; ++i) {
        U[i] = P[2 * i];
        V[i] = P[2 * i + 1];
    }
}

void combine_results_FFT(std::vector<complex>& P, const std::vector<complex>& U, const std::vector<complex>& V, size_t N, size_t start, size_t end) {
    for (size_t i = start; i < end; ++i) {
        complex t = std::polar(1.0, -2 * M_PI * i / N);
        P[i] = U[i] + t * V[i];
        P[i + N / 2] = U[i] - t * V[i];
    }
}

void FFT_Radix2_Parallel(std::vector<complex>& P, size_t num_threads) {
    size_t N = P.size();
    if (!is_power_of_two(N)) {
        std::cout << "The input size " << N << " is not a power of 2! Resizing input" << std::endl;
        N = next_power_of_two(N);
        P.resize(N);
    }

    if (num_threads == 1 || N == 1) {
        P = Radix2FFT(P);
        return;
    }

    std::vector<complex> U(N / 2), V(N / 2);

    size_t split_threads = std::min(num_threads, N / 2);
    size_t chunk_size = (N / 2 + split_threads - 1) / split_threads;
    std::vector<std::thread> split_threads_list;

    for (size_t t = 0; t < split_threads; ++t) {
        size_t start = t * chunk_size;
        size_t end = std::min(start + chunk_size, N / 2);
        split_threads_list.emplace_back(split_data_FFT, std::cref(P), std::ref(U), std::ref(V), start, end);
    }

    for (auto& th : split_threads_list) {
        th.join();
    }

    std::thread t1(FFT_Radix2_Parallel, std::ref(U), num_threads / 2);
    std::thread t2(FFT_Radix2_Parallel, std::ref(V), num_threads / 2);
    t1.join();
    t2.join();

    split_threads_list.clear();
    for (size_t t = 0; t < split_threads; ++t) {
        size_t start = t * chunk_size;
        size_t end = std::min(start + chunk_size, N / 2);
        split_threads_list.emplace_back(combine_results_FFT, std::ref(P), std::cref(U), std::cref(V), N, start, end);
    }

    for (auto& th : split_threads_list) {
        th.join();
    }
}



/////////////////////////////// IFFT //////////////////////////////////////////

void split_data_IFFT(const std::vector<complex>& P_star, std::vector<complex>& U, std::vector<complex>& V, size_t start, size_t end) {
    for (size_t i = start; i < end; ++i) {
        U[i] = P_star[2 * i];
        V[i] = P_star[2 * i + 1];
    }
}


void combine_results_IFFT(std::vector<complex>& P_star, const std::vector<complex>& U, const std::vector<complex>& V, size_t N, size_t start, size_t end) {
    for (size_t i = start; i < end; ++i) {
        complex t = std::polar(1.0, -2 * M_PI * i / N);
        t = std::conj(t);
        P_star[i] = (U[i] + t * V[i]) / 2.0;
        P_star[i + N / 2] = (U[i] - t * V[i]) / 2.0;
    }
}

void Inverse_Radix2_Parallel(std::vector<complex>& P_star, size_t num_threads) {

    size_t N = P_star.size();
    if (!is_power_of_two(N)) {
        std::cout << "The input size " << N << " is not a power of 2! Resizing input" << std::endl;
        N = next_power_of_two(N);
        P_star.resize(N);
    }

    if (num_threads == 1 || N == 1) {
        P_star = InverseRadix2FFT(P_star);
        return;
    }

    std::vector<complex> U(N / 2), V(N / 2);

    size_t split_threads = std::min(num_threads, N / 2);
    size_t chunk_size = (N / 2 + split_threads - 1) / split_threads;
    std::vector<std::thread> split_threads_list;

    for (size_t t = 0; t < split_threads; ++t) {
        size_t start = t * chunk_size;
        size_t end = std::min(start + chunk_size, N / 2);
        split_threads_list.emplace_back(split_data_IFFT, std::cref(P_star), std::ref(U), std::ref(V), start, end);
    }

    for (auto& th : split_threads_list) {
        th.join();
    }

    std::thread t1(Inverse_Radix2_Parallel, std::ref(U), num_threads / 2);
    std::thread t2(Inverse_Radix2_Parallel, std::ref(V), num_threads / 2);
    t1.join();
    t2.join();

    split_threads_list.clear();
    for (size_t t = 0; t < split_threads; ++t) {
        size_t start = t * chunk_size;
        size_t end = std::min(start + chunk_size, N / 2);
        split_threads_list.emplace_back(combine_results_IFFT, std::ref(P_star), std::cref(U), std::cref(V), N, start, end);
    }

    for (auto& th : split_threads_list) {
        th.join();
    }
}