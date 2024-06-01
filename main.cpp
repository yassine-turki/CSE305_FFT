#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <fstream>
#include <string>
#include <thread>
#include <chrono>
#include <sstream>
#include <stdexcept>

#define M_PI 3.14159265358979323846

typedef std::complex<double> complex;



bool is_power_of_two(int n){
    return (n & (n - 1)) == 0;
}

int next_power_of_two(int n){
    int m = std::ceil(std::log2(n));
    return std::pow(2, m);
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

std::vector<complex> get_weather_data(std::string file_path) {
    std::ifstream file(file_path); // Open the CSV file
    if (!file.is_open()) {
        std::cerr << "Failed to open the file!" << std::endl;

    }

    std::vector<complex> weather_vector; // Vector to store the first column values
    std::string line;
    // Skip the first row (header)
    if (!std::getline(file, line)) {
        std::cerr << "Failed to read the header line!" << std::endl;

    }

    // Read the file line by line
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string date;
        std::string value;

        // Extract the first column value
        if (std::getline(ss, date, ',')) {
            if(std::getline(ss, value, ',')) {
                double val = std::stod(value);
                weather_vector.push_back(val);
            }
        }
    }

    file.close(); // Close the file
    weather_vector.resize(next_power_of_two(weather_vector.size()));
    return weather_vector;

}

// Polynomial operations

complex evaluate_polynomial_trivial(std::vector<complex> P, complex x){
    int N = P.size();
    complex X = 1.;
    complex y;
    for (int j = 0; j < N; j++){
        y += P[j] * X;
        X *= x;
    }
    return y;
}

std::vector<complex> add_polynomials(std::vector<complex> P, std::vector<complex> Q){
    int N = std::max(P.size(), Q.size());
    P.resize(N);
    Q.resize(N);
    std::vector<complex> R(N);
    for (int j = 0; j < N; j++){
        R[j] = P[j] + Q[j];
    }
    return R;
}

std::vector<complex> multiply_polynomials_trivial(std::vector<complex> P, std::vector<complex> Q){
    int N = P.size();
    int M = Q.size();
    std::vector<complex> R(M + N, 0.);
    for (int j = 0; j < N; j++){
        for (int k = 0; k < M; k++){
            R[j + k] += P[j] * Q[k];
        }
    }
    return R;
}


complex horner_evaluate(std::vector<complex> P, complex x){
    int N = P.size();
    complex y = P[N - 1];
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

std::vector<complex> generate_synthetic_data(int num_points) {
    std::vector<complex> data(num_points);
    double frequency = 1.0;
    double amplitude = 1.5;
    for (int i = 0; i < num_points; ++i) {
        data[i] = amplitude * std::sin(2 * M_PI * frequency * i / num_points);
    }
    return data;
}



int main() {

    bool test_alex = 0;
    std::string file_path;
    if(test_alex){
//         file_path = "data/natural_gas_co2_emissions_for_electric_power_sector.csv";
         file_path = "data/Paris_data.csv";
        std::vector<complex> weather_data = get_weather_data(file_path);

    }
    else {

        // Test with dummy data
//        std::vector<complex> weather_data = {1, 2, 3, 4};
//        std::vector<complex> weather_data_parallel = {1, 2, 3, 4};
//        std::vector<complex> result_seq = {1, 2, 3, 4};
//        std::vector<complex> result_parallel = {1, 2, 3, 4};

        // Test with data from file

//        file_path = "../data/natural_gas_co2_emissions_for_electric_power_sector.csv";
           file_path = "../data/Paris_data.csv";


        std::vector<complex> weather_data = get_weather_data(file_path);
        std::vector<complex> weather_data_parallel =  get_weather_data(file_path);
        std::vector<complex> result_seq = get_weather_data(file_path);
        std::vector<complex> result_parallel =  get_weather_data(file_path);


//        std::vector<complex> weather_data = generate_synthetic_data(2048);

//


//        std::cout << "Size of the input vector: " << weather_data.size() << std::endl;

        std::chrono :: time_point<std::chrono::system_clock> start_seq, end_seq, start_parallel, end_parallel;
        start_seq = std::chrono::system_clock::now();
        result_seq = Radix2FFT(weather_data);

        // Lines to check IFFT
        result_seq = InverseRadix2FFT(result_seq);

        end_seq = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds_seq = end_seq - start_seq;
        std::cout << "Time taken for serial FFT: " << elapsed_seconds_seq.count() << "s\n";


        int num_threads = 10;
        start_parallel = std::chrono::system_clock::now();
        FFT_Radix2_Parallel(result_parallel, num_threads);

        // Lines to check IFFT
        Inverse_Radix2_Parallel(result_parallel, num_threads);
        end_parallel = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds_parallel = end_parallel - start_parallel;

        std::cout << "Time taken for parallel FFT with " << num_threads << " threads: " << elapsed_seconds_parallel.count() << "s\n";

        // Comparison of Sequential and parallel
        for (int i = 0; i < weather_data.size(); i++) {
            if (std::abs(result_seq[i] - result_parallel[i]) > 1e-6){
                std::cout << "Error at index " << i << std::endl;
                std::cout << "Sequential: " << result_seq[i] << " Parallel: " << result_parallel[i] << std::endl;
            }
        }
        // Comparison of sequential with actual data
        for (int i = 0; i < weather_data.size(); i++) {
            if (std::abs(result_seq[i] - weather_data[i]) > 1e-6){
                std::cout << "Error at index " << i << std::endl;
                std::cout << "Sequential: " << result_seq[i] << " Actual: " << weather_data[i] << std::endl;
            }
        }

        // Comparison of parallel with actual data
        for (int i = 0; i < weather_data.size(); i++) {
            if (std::abs(result_parallel[i] - weather_data_parallel[i]) > 1e-6){
                std::cout << "Error at index " << i << std::endl;
                std::cout << "Parallel: " << result_parallel[i] << " Actual: " << weather_data_parallel[i] << std::endl;
            }
        }


    }

    return 0;
}