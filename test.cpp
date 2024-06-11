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
#include "benchmark.cpp"

#define M_PI 3.14159265358979323846

typedef std::complex<double> complex;



/*------------------------------------Testing functions ------------------------------------*/



std::vector<double> test_naive_dft(std::vector<complex> data, int num_threads) {
    std::vector<double> times;
    std::vector<complex> result_seq, result_parallel;
    auto start_seq = std::chrono::system_clock::now();

    result_seq = Naive_DFT(data);
    result_seq = Naive_IDFT(result_seq);

    auto end_seq = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds_seq = end_seq - start_seq;
    std::cout << "Time taken for Naive serial DFT: " << elapsed_seconds_seq.count() << "s\n";
    times.push_back(elapsed_seconds_seq.count());

    auto start_parallel = std::chrono::system_clock::now();

    result_parallel = Naive_DFT_Parallel(data, num_threads);
    result_parallel = Naive_IDFT_Parallel(result_parallel, num_threads);

    auto end_parallel = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds_parallel = end_parallel - start_parallel;
    std::cout << "Time taken for Naive parallel DFT with " << num_threads << " threads: " << elapsed_seconds_parallel.count() << "s\n";
    times.push_back(elapsed_seconds_parallel.count());

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

        return times;
}

std::vector<double> test_radix2(std::vector<complex> data, int num_threads) {
    std::vector<double> times;
    std::vector<complex> result_seq;
    std::vector<complex> result_parallel = data;
    auto start_seq = std::chrono::system_clock::now();

    result_seq = Radix2FFT(data);
    result_seq = InverseRadix2FFT(result_seq);

    auto end_seq = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds_seq = end_seq - start_seq;
    std::cout << "Time taken for serial Radix2: " << elapsed_seconds_seq.count() << "s\n";
    times.push_back(elapsed_seconds_seq.count());

    auto start_parallel = std::chrono::system_clock::now();

    FFT_Radix2_Parallel(result_parallel, num_threads);
    Inverse_Radix2_Parallel(result_parallel, num_threads);

    auto end_parallel = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds_parallel = end_parallel - start_parallel;
    std::cout << "Time taken for parallel Radix2 with " << num_threads << " threads: " << elapsed_seconds_parallel.count() << "s\n";
    times.push_back(elapsed_seconds_parallel.count());

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

        return times;
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

void test_prime_decomposition(int n) {
    std::vector<int> prime_factors = prime_decomposition(n);
    std::cout << "Testing n: " << n << "\n";
    for (const auto& c : prime_factors) {
        std::cout << c <<" ";
    }
    std::cout<<"Sequential"<<std::endl;
    int num_threads = 5;
    std::vector<int> prime_factors_2 = prime_decomposition_parallel(n, num_threads);
    for (const auto& c : prime_factors_2) {
        std::cout << c <<" ";
    }
    std::cout<<"Parallel"<<std::endl;
}

std::vector<double> test_convolutions(std::pair<std::vector<int>, std::vector<int>> &data, std::unordered_map<int, std::pair<int, int>> roots, int num_threads){
    std::vector<double> times;
    std::vector<int> poly_1 = data.first;
    std::vector<int> poly_2 = data.second;
    int n = poly_1.size();
    int p = prime_ntt_find(n);

    // Naive ntt
    // auto start_naive = std::chrono::high_resolution_clock::now();
    // std::cout<<"Enter 1"<<std::endl;
    // std::vector<int> conv_naive = convolution_ntt(poly_1, poly_2, p, roots);
    // auto end_naive = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double> duration_naive = end_naive - start_naive;
    // times.push_back(duration_naive.count());

    // Fast ntt
    // auto start_fast = std::chrono::high_resolution_clock::now();
    // std::cout<<"Enter 2"<<std::endl;
    // std::vector<int> conv_fast = convolution_fast(poly_1, poly_2, p, roots);
    // auto end_fast = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double> duration_fast = end_fast - start_fast;
    // times.push_back(duration_fast.count());

    // Naive ntt parallel
    // auto start_naive_parallel = std::chrono::high_resolution_clock::now();
    // std::cout<<"Enter 3"<<std::endl;
    // std::vector<int> conv_naive_parallel = convolution_ntt_parallel(poly_1, poly_2, p, roots, num_threads);
    // auto end_naive_parallel = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double> duration_naive_parallel = end_naive_parallel - start_naive_parallel;
    // times.push_back(duration_naive_parallel.count());

    // Fast ntt parallel
    auto start_fast_parallel = std::chrono::high_resolution_clock::now();
    std::cout<<"Enter 4"<<std::endl;
    std::vector<int> conv_fast_parallel = convolution_fast_parallel(poly_1, poly_2, p, roots, num_threads);
    auto end_fast_parallel = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration_fast_parallel = end_fast_parallel - start_fast_parallel;
    times.push_back(duration_fast_parallel.count());   
    return times;    
}

// Function to remove the extension from a file path

std::string removeExtension(std::string& imagePath) {
    size_t lastDot = imagePath.find_last_of(".");
    if (lastDot == std::string::npos) {
        return imagePath;
    }
    return imagePath.substr(0, lastDot);
}

void test_compression_weather_data(std::vector<complex> data, std::vector<int> percentages,size_t num_threads){
    for(auto i : percentages){
        int percentage_keep = i;
        std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
        FFT_Radix2_Parallel(data, num_threads);
        double threshold = get_threshold_value(data, percentage_keep/100.);
        compress(data, threshold, num_threads);
        Inverse_Radix2_Parallel(data, num_threads);
        std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end - start;
        std::string percentage_keep_str = std::to_string(percentage_keep);
        std::string outfile_path = "../data/weather_data_keep" + percentage_keep_str + ".csv";
        write_data(data, outfile_path);
        std::cout << "Time taken for compression with " << percentage_keep << "% kept: " << duration.count() << "s\n";
    }
}

void test_compression_image(std::string image_path, std::vector<int> percentages, size_t num_threads){
    int w,h;
    for(auto i : percentages){
        int percentage_keep = i;
        std::vector<complex> *img_data = read_img(image_path, w, h);
        std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now();
        compute_fft_img(*img_data, w, h, num_threads);
        double threshold = get_threshold_value(*img_data, percentage_keep/100.);
        compress(*img_data, threshold, num_threads);
        compute_ifft_img(*img_data, w, h, num_threads);
        std::string percentage_keep_str = std::to_string(percentage_keep);
        std::string outfile_path = removeExtension(image_path) + "_keep" + percentage_keep_str + ".png";
        write_img(*img_data, w, h, outfile_path);
        std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration = end - start;
        std::cout << "Time taken for compression with " << percentage_keep << "% kept: " << duration.count() << "s\n";
    }

}