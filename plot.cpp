#include <iostream>
#include <vector>
#include <complex>
#include "test.cpp"


std::vector<double> benchmark_naive_dft(std::vector<complex> data, std::vector<int> num_threads, int num_executions) {
    std::vector<double> times;
    std::vector<complex> result_seq;
    std::vector<complex> result_parallel = data;

    for (int i = 0; i < num_executions; ++i) {
        double seq_time_sum = 0.0;
        std::vector<double> parallel_time_sum(num_threads.size(), 0.0);

        auto start_seq = std::chrono::system_clock::now();
        result_seq = Naive_DFT(data);
        result_seq = Naive_IDFT(result_seq);
        auto end_seq = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds_seq = end_seq - start_seq;
        seq_time_sum += elapsed_seconds_seq.count();

        for (int j = 0; j < num_threads.size(); ++j) {
            auto start_parallel = std::chrono::system_clock::now();
            Naive_DFT_Parallel(result_parallel, num_threads[j]);
            Naive_IDFT_Parallel(result_parallel, num_threads[j]);
            auto end_parallel = std::chrono::system_clock::now();
            std::chrono::duration<double> elapsed_seconds_parallel = end_parallel - start_parallel;
            parallel_time_sum[j] += elapsed_seconds_parallel.count();
        }

        times.push_back(seq_time_sum);
        for (double time : parallel_time_sum) {
            times.push_back(time);
        }
    }

    for (int i = 0; i < times.size(); ++i) {
        times[i] /= num_executions;
    }

    return times;
}

std::vector<double> benchmark_radix2(std::vector<complex> data, std::vector<int> num_threads, int num_executions) {
    std::vector<double> times;
    std::vector<complex> result_seq;
    std::vector<complex> result_parallel = data;

    for (int i = 0; i < num_executions; ++i) {
        double seq_time_sum = 0.0;
        std::vector<double> parallel_time_sum(num_threads.size(), 0.0);

        auto start_seq = std::chrono::system_clock::now();
        result_seq = Radix2FFT(data);
        result_seq = InverseRadix2FFT(result_seq);
        auto end_seq = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds_seq = end_seq - start_seq;
        seq_time_sum += elapsed_seconds_seq.count();

        for (int j = 0; j < num_threads.size(); ++j) {
            auto start_parallel = std::chrono::system_clock::now();
            FFT_Radix2_Parallel(result_parallel, num_threads[j]);
            Inverse_Radix2_Parallel(result_parallel, num_threads[j]);
            auto end_parallel = std::chrono::system_clock::now();
            std::chrono::duration<double> elapsed_seconds_parallel = end_parallel - start_parallel;
            parallel_time_sum[j] += elapsed_seconds_parallel.count();
        }

        times.push_back(seq_time_sum);
        for (double time : parallel_time_sum) {
            times.push_back(time);
        }
    }

    for (int i = 0; i < times.size(); ++i) {
        times[i] /= num_executions;
    }

    return times;
}


int benchmark(std::vector<complex> data, int num_executions = 10){

    // Every time vector contains a pair of values: time for sequential and parallel
    std::vector<int> num_threads;
    for(int i = 2; i <= 25; i+=2) {
        num_threads.push_back(i);
    }
    std::vector<double> times_naive_dft = benchmark_naive_dft(data, num_threads, num_executions);
    std::cout << "Naive DFT benchmark done." << std::endl;
    std::vector<double> times_radix2 = benchmark_radix2(data, num_threads, num_executions);
    std::cout << "Radix2 benchmark done." << std::endl;

    std::ofstream file("../execution_times.csv");
    if (!file.is_open()) {
        std::cerr << "Failed to open file." << std::endl;
        return 1;
    }

    // Write header
    file << "Algorithm";
    for (int num_thread : num_threads) {
        file << "," << num_thread;
    }
    file << "\n";

    file << "Naive DFT";

    double sequential_time_naive_dft = times_naive_dft[0];
    file << "," << sequential_time_naive_dft;
    for (int i=1; i < num_threads.size(); i++) {
        double time = times_naive_dft[i];
        file << "," << time;
    }
    file << "\n";

    file << "Radix2";
    double sequential_time_radix2 = times_radix2[0];
    file << "," << sequential_time_radix2;
    for (int i=1; i < num_threads.size(); i++) {
        double time = times_radix2[i];
        file << "," << time;
    }
    file << "\n";

    file.close();
    std::cout << "CSV file generated successfully." << std::endl;
    return 0;


}
