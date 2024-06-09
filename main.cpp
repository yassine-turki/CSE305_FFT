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
#include "compress.cpp"


#define M_PI 3.14159265358979323846

typedef std::complex<double> complex;



int main() {

    bool test_mult = 0;
    bool test_fft = 0;
    bool test_compression = 1;
    std::string file_path;
    std::vector<complex> data_vect;
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

        // Test compression image


//        std::string image_path =  "../philip.png";
//        int w,h;
//        std::vector<complex> *img_data = read_img(image_path, w, h);
//        std::vector<complex> img_data_copy = *img_data;
//        int num_threads = 10;
//        parallel_fft_image(*img_data, w, h, num_threads);
//       double threshold = get_threshold_value(*img_data, 0.8);
//
//        double threshold = get_threshold_value(*img_data, 0.4);
//        compress(*img_data, threshold, num_threads);
//        parallel_inverse_fft_image(*img_data, w, h, num_threads);
//        write_img(*img_data, w, h, "../test.png");

// Test compression weather_data

      std::vector<complex> weather_data = get_weather_data("../data/weather_data_2024.csv");
      int num_threads = 10;
      FFT_Radix2_Parallel(weather_data, num_threads);
      for(int i = 40; i>0; i-=39){
            int percentage_keep = i;
            double threshold = get_threshold_value(weather_data, percentage_keep / 100.);
            compress(weather_data, threshold, num_threads);
            Inverse_Radix2_Parallel(weather_data, num_threads);
            std::string percentage_keep_str = std::to_string(percentage_keep);
            std::string compressed_data_path = "../data/weather_data_2024_" + percentage_keep_str + ".csv";
            write_data(weather_data, compressed_data_path);
        }




// Write data synthetic
//
//       std::vector<complex> synthetic_data = generate_synthetic_data(10E4);
//       write_data(synthetic_data, "../data/synthetic_data10E4.csv");
//       int num_threads = 10;
//       FFT_Radix2_Parallel(synthetic_data, num_threads);
//       int percentage_keep = 40;
//       double threshold = get_threshold_value(synthetic_data, percentage_keep/100.);
//       compress(synthetic_data, threshold, num_threads);
//       Inverse_Radix2_Parallel(synthetic_data, num_threads);
//       std::string percentage_keep_str = std::to_string(percentage_keep);
//       std::string compressed_data_path = "../data/synthetic_data10E4_" + percentage_keep_str + ".csv";
//       write_data(synthetic_data, compressed_data_path);
//        percentage_keep = 1;
//        threshold = get_threshold_value(synthetic_data, percentage_keep/100.);
//        compress(synthetic_data, threshold, num_threads);
//        Inverse_Radix2_Parallel(synthetic_data, num_threads);
//        percentage_keep_str = std::to_string(percentage_keep);
//        compressed_data_path = "../data/synthetic_data10E4_" + percentage_keep_str + ".csv";
//        write_data(synthetic_data, compressed_data_path);


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