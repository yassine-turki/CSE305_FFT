#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <fstream>
#include <string>
#include "fft_multiplication.cpp"

#define M_PI 3.14159265358979323846

typedef std::complex<double> complex;


/*------------------------------------Data functions ------------------------------------*/

std::vector<complex> generate_dummy_data() {
    return {1, 2, 3, 4};
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


std::vector<complex> generate_synthetic_data(int num_points) {
    std::vector<complex> data(num_points);
    double frequency = 5;
    double amplitude = 10;
    for (int i = 0; i < num_points; ++i) {
        data[i] = amplitude * std::sin(2 * M_PI * frequency * i / num_points);
    }
    return data;
}




void write_data(std::vector<complex> &data, std::string file_path) {
std::ofstream file(file_path);


    if (!file.is_open()) {
        std::cerr << "Failed to open the file!" << std::endl;
    }
    // header
    file << "Value\n";

    for (int i = 0; i < data.size(); i++) {
        file << data[i].real() << "\n";
    }
    file.close();
}