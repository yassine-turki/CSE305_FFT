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
#include <unordered_map>
#include <random>
#include <mutex>
#include <tuple>
#include "test.cpp"


#define M_PI 3.14159265358979323846

typedef std::complex<double> complex;

int main() {

    bool test_mult = 1; // Test polynomial multiplication
    bool test_fft = 0; // Test FFT Algorithms
    bool test_compression = 0; // Test Compression
    std::string file_path;
    std::vector<complex> data_vect;
    if(test_mult){
        // int p = 7681;
        std::vector<int> intVector = {1, 2, 3, 4};
        std::vector<int> intVector2 = {5, 6, 7, 8};
        int num_threads = 13;
        std::unordered_map<int, std::pair<int, int>> roots = precomputed_roots();
        std::vector<int> inverse_ntt = {2489, 7489, 6478, 6607};
        std::pair<std::vector<int>, std::vector<int>> data = generate_synthetic_data_ntt(16384);
        // int p = prime_ntt_find(4096);
       
        // std::vector<int> conv_naive = convolution_ntt(data.first, data.second, p, roots);
        // for (const auto& term: conv_naive){
        //     std::cout<<term<<" ";
        // }
        // std::cout<<"Naive"<<std::endl;

        // std::vector<int> conv_fast = convolution_fast(data.first, data.second, p, roots);
        // for (const auto& term: conv_fast){
        //     std::cout<<term<<" ";
        // }
        // std::cout<<"Fast"<<std::endl;

        // std::vector<int> conv_naive_parallel = convolution_ntt_parallel(data.first, data.second, p, roots, num_threads);
        // for (const auto& term: conv_naive_parallel){
        //     std::cout<<term<<" ";
        // }
        // std::cout<<"Naive parallel"<<std::endl;

        // std::vector<int> conv_fast_parallel = convolution_fast_parallel(data.first, data.second, p, roots, num_threads);
        // for (const auto& term: conv_fast_parallel){
        //     std::cout<<term<<" ";
        // }
        // std::cout<<"Fast parallel"<<std::endl;

        std::vector<double> times = test_convolutions(data, roots, num_threads);
        // std::cout<<"Time taken for naive convolution is "<<times[0]<<std::endl;
        // std::cout<<"Time taken for fast convolution is "<<times[0]<<std::endl;
        // std::cout<<"Time taken for naive parallel convolution with "<<num_threads<<" threads is "<<times[1]<<std::endl;
        std::cout<<"Time taken for fast parallel convolution with "<<num_threads<<" threads is "<<times[0]<<std::endl;

    }
    else if(test_fft) {

        bool use_dummy_data = 0; // Test it on {1, 2, 3, 4}
        bool use_synthetic_data = 1; // Test it on synthetic data (sine)
        bool use_weather_data = 0; // Test it on weather data from file
        int num_threads = 10;

        if (use_dummy_data) {
            data_vect = {1, 2, 3, 4};
        } else if (use_synthetic_data) {
            data_vect = generate_synthetic_data(10E4);
        } else if (use_weather_data) {
            file_path = "../data/natural_gas_co2_emissions_for_electric_power_sector.csv";
            //file_path = "../data/Paris_data.csv";
            data_vect = get_weather_data(file_path);
        }

//        test_naive_dft(data_vect, num_threads);
//        test_radix2(data_vect, num_threads);
        int num_executions = 5;
        benchmark(data_vect, num_executions);
    }
    else if(test_compression){
        bool test_image = 0;
        bool test_weather_data = 1;
        int num_threads = 20;
        std::vector<int> percentages = {70, 50, 40, 30, 20, 10, 1};

        if(test_image){
        // Test compression image

        // MONA LISA
        std::string image_path =  "../images/mona_lisa/monalisa.jpg";
        test_compression_image(image_path, percentages, num_threads);

        // HORSE IMAGE
       image_path =  "../images/horse/horse.jpg";
       test_compression_image(image_path, percentages, num_threads);


        // Test compression weather_data

        std::vector<complex> weather_data = get_weather_data("../data/weather_data_2024.csv");
        test_compression_weather_data(weather_data, percentages, num_threads);

    }
        if(test_weather_data){
            // Test compression weather_data

            std::vector<complex> weather_data = get_weather_data("../data/weather_data_2024.csv");
            test_compression_weather_data(weather_data, percentages, num_threads);
        }

    else{

        // Test is_prime

//        size_t num_threads = 5;
//        std::vector<int> test_values = {29, 97, 100, 101, 999983, 1000003, 36, 144, 185, 17};
//        test_is_prime(test_values, num_threads);

// Test find_generator

//        size_t num_threads = 5;
//        int n = 16;
//        int p = 7681;
//        test_find_generator(n, p, num_threads);



    }

    return 0;
}
}