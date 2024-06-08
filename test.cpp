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
#include "data.cpp"

#define M_PI 3.14159265358979323846

typedef std::complex<double> complex;



/*------------------------------------Testing functions ------------------------------------*/



void test_naive_dft(std::vector<complex> data, int num_threads) {
    std::vector<complex> result_seq, result_parallel;
    auto start_seq = std::chrono::system_clock::now();

    result_seq = Naive_DFT(data);
    result_seq = Naive_IDFT(result_seq);

    auto end_seq = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds_seq = end_seq - start_seq;
    std::cout << "Time taken for Naive serial DFT: " << elapsed_seconds_seq.count() << "s\n";

    auto start_parallel = std::chrono::system_clock::now();

    result_parallel = Naive_DFT_Parallel(data, num_threads);
    result_parallel = Naive_IDFT_Parallel(result_parallel, num_threads);

    auto end_parallel = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds_parallel = end_parallel - start_parallel;
    std::cout << "Time taken for Naive parallel DFT with " << num_threads << " threads: " << elapsed_seconds_parallel.count() << "s\n";

    // Comparison of Sequential with data
    bool error = false;
//    std::cout << "Comparing Naive DFT sequential with data" << std::endl;
    for (size_t i = 0; i < data.size(); i++) {
        if (std::abs(result_seq[i] - data[i]) > 1e-6) {
            std::cout << "Error at index " << i << std::endl;
            std::cout << "Data: " << data[i] << " Sequential: " << result_seq[i] << std::endl;
            error = true;
        }
    }
//    if(!error) {
//        std::cout << "No errors found!" << std::endl;
//    }

    error = false;
//    std::cout << "Comparing Naive DFT parallel with data" << std::endl;
    for (size_t i = 0; i < data.size(); i++) {
        if (std::abs(result_parallel[i] - data[i]) > 1e-6) {
            std::cout << "Error at index " << i << std::endl;
            std::cout << "Data: " << data[i] << " Parallel: " << result_parallel[i] << std::endl;
            error = true;
        }
    }
//    if(!error) {
//        std::cout << "No errors found!" << std::endl;
//    }
}

void test_radix2(std::vector<complex> data, int num_threads) {
    std::vector<complex> result_seq;
    std::vector<complex> result_parallel = data;
    auto start_seq = std::chrono::system_clock::now();

    result_seq = Radix2FFT(data);
    result_seq = InverseRadix2FFT(result_seq);

    auto end_seq = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds_seq = end_seq - start_seq;
    std::cout << "Time taken for serial Radix2: " << elapsed_seconds_seq.count() << "s\n";

    auto start_parallel = std::chrono::system_clock::now();

    FFT_Radix2_Parallel(result_parallel, num_threads);
    Inverse_Radix2_Parallel(result_parallel, num_threads);

    auto end_parallel = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds_parallel = end_parallel - start_parallel;
    std::cout << "Time taken for parallel Radix2 with " << num_threads << " threads: " << elapsed_seconds_parallel.count() << "s\n";

    // Comparison of Sequential with data
    bool error = false;
//    std::cout << "Comparing Radix2 sequential with data" << std::endl;
    for (size_t i = 0; i < data.size(); i++) {
        if (std::abs(result_seq[i] - data[i]) > 1e-6) {
            std::cout << "Error at index " << i << std::endl;
            std::cout << "Data: " << data[i] << " Sequential: " << result_seq[i] << std::endl;
            error = true;
        }
    }
//    if(!error) {
//        std::cout << "No errors found!" << std::endl;
//    }

    error = false;
//    std::cout << "Comparing Radix2 parallel with data" << std::endl;
    for (size_t i = 0; i < data.size(); i++) {
        if (std::abs(result_parallel[i] - data[i]) > 1e-6) {
            std::cout << "Error at index " << i << std::endl;
            std::cout << "Data: " << data[i] << " Parallel: " << result_parallel[i] << std::endl;
            error = true;
        }
    }
//    if(!error) {
//        std::cout << "No errors found!" << std::endl;
//    }
}


void test_is_prime(const std::vector<int>& test_values, size_t num_threads) {
    for (int p : test_values) {
        auto start = std::chrono::high_resolution_clock::now();
        bool result_single = is_prime(p);
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration_single = end - start;

        start = std::chrono::high_resolution_clock::now();
        bool result_parallel = is_prime_parallel(p, num_threads);
        end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration_parallel = end - start;

        std::cout << "Testing number: " << p << "\n";
        std::cout << "Single-threaded result: " << result_single << ", Time: " << duration_single.count() << " seconds\n";
        std::cout << "Multi-threaded result with " << num_threads << " threads: "<<result_parallel << ", Time: " << duration_parallel.count() << " seconds\n";
        std::cout << "--------------------------------------------------\n";
    }
}