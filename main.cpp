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
#include "test.cpp"


#define M_PI 3.14159265358979323846

typedef std::complex<double> complex;

int main() {

    bool test_mult = 0;
    bool test_fft = 0;
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

        std::vector<int> int_arr = convolution_ntt(intVector, intVector2);
        for (const auto& c : int_arr) {
            std::cout << c <<" ";
        }
        std::cout<<"NTT"<<std::endl;
        std::vector<int> int_arr_2 = convolution_ntt(intVector, intVector2);
        for (const auto& c : int_arr_2) {
            std::cout << c <<" ";
        }
        // std::cout<<"fast NTT"<<std::endl;
        // int num_threads = 10;
        // std::vector<int> int_sum_3 = ntt_parallel(intVector, p, num_threads);
        // for (const auto& c : int_sum_3) {
        //     std::cout << c <<" ";
        // }
        // std::cout<<"parallel NTT with "<< num_threads<<" threads"<<std::endl;

        // std::vector<int> int_sum_4 = intVector;
        // fast_ntt_parallel(int_sum_4, p, num_threads);
        // for (const auto& c : int_sum_4) {
        //     std::cout << c <<" ";
        // }
        // std::cout<<"parallel fast NTT with "<< num_threads<<" threads"<<std::endl;


    }
    else if(test_fft) {

        bool use_dummy_data = 0; // Test it on {1, 2, 3, 4}
        bool use_synthetic_data = 0; // Test it on synthetic data (sine)
        bool use_weather_data = 1; // Test it on weather data from file
        int num_threads = 10;

        if (use_dummy_data) {
            weather_data = {1, 2, 3, 4};
        } else if (use_synthetic_data) {
            weather_data = generate_synthetic_data(10E4);
        } else if (use_weather_data) {
            file_path = "../data/natural_gas_co2_emissions_for_electric_power_sector.csv";
            //file_path = "../data/Paris_data.csv";
            weather_data = get_weather_data(file_path);
        }

        test_naive_dft(weather_data, num_threads);
        test_radix2(weather_data, num_threads);
    }
    else{

        // Test is_prime

//        size_t num_threads = 5;
//        std::vector<int> test_values = {29, 97, 100, 101, 999983, 1000003, 36, 144, 185, 17};
//        test_is_prime(test_values, num_threads);

// Test find_generator

        size_t num_threads = 5;
        int n = 16;
        int p = 7681;
        test_find_generator(n, p, num_threads);
    }


    return 0;
}