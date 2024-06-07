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
#include "fft_multiplication.cpp"


#define M_PI 3.14159265358979323846

typedef std::complex<double> complex;


/*------------------------------------Data functions ------------------------------------*/

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


std::vector<complex> generate_synthetic_data(int num_points) {
    std::vector<complex> data(num_points);
    double frequency = 1.0;
    double amplitude = 1.5;
    for (int i = 0; i < num_points; ++i) {
        data[i] = amplitude * std::sin(2 * M_PI * frequency * i / num_points);
    }
    return data;
}


/*------------------------------------Testing functions ------------------------------------*/

std::vector<complex> generate_dummy_data() {
    return {1, 2, 3, 4};
}

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

int main() {

    bool test_mult = 1;
    std::string file_path;
    std::vector<complex> weather_data;
    if(test_mult){

        std::vector<int> intVector = {1, 2, 3, 4};
        std::vector<int> intVector2 = {5, 6, 7, 8};
        int p = 7681;
        std::vector<int> inverse_ntt = {2489, 7489, 6478, 6607};

        std::vector<complex> complexVector = {
                {1.0, 2.0},
                {3.0, 4.0},
                {5.0, 6.0},
                {7.0, 8.0}
        };
        std::vector<complex> complexVector2 = {
                {1.0, 2.0},
                {3.0, 4.0},
                {5.0, 6.0},
                {7.0, 8.0}
        };

        std::vector<complex> complexVector3 = {
                {1.0, 0.0},
                {2.0, 0.0},
                {3.0, 0.0},
                {4.0, 0.0}
        };
        std::vector<complex> complexVector4 = {
                {1.0, 0.0},
                {2.0, 0.0},
                {3.0, 0.0},
                {4.0, 0.0}
        };

        std::vector<int> int_sum = ntt(intVector, p);
        for (const auto& c : int_sum) {
            std::cout << c <<" ";
        }
        std::cout<<"FFT"<<std::endl;
        std::vector<int> int_sum_2 = fast_ntt(intVector, p);
        for (const auto& c : int_sum_2) {
            std::cout << c <<" ";
        }

    }
    else {

        bool use_dummy_data = 0; // Test it on {1, 2, 3, 4}
        bool use_synthetic_data = 0; // Test it on synthetic data (sine)
        bool use_weather_data = 1; // Test it on weather data from file
        int num_threads = 10;

        if (use_dummy_data) {
            weather_data = {1, 2, 3, 4};
        } else if (use_synthetic_data) {
            weather_data = generate_synthetic_data(10E4);
        } else if (use_weather_data) {
            file_path = "data/natural_gas_co2_emissions_for_electric_power_sector.csv";
            //file_path = "../data/Paris_data.csv";
            weather_data = get_weather_data(file_path);
        }

        test_naive_dft(weather_data, num_threads);
        test_radix2(weather_data, num_threads);
    }


    return 0;
}