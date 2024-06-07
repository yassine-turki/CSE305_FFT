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

        std::vector<int> int_sum = intt(inverse_ntt, p);
        for (const auto& c : int_sum) {
            std::cout << c <<" ";
        }
        std::cout<<"FFT"<<std::endl;
        std::vector<int> int_sum_2 = fast_intt(inverse_ntt, p);
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
            file_path = "../data/natural_gas_co2_emissions_for_electric_power_sector.csv";
            //file_path = "../data/Paris_data.csv";
            weather_data = get_weather_data(file_path);
        }

        test_naive_dft(weather_data, num_threads);
        test_radix2(weather_data, num_threads);
    }


    return 0;
}