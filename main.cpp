#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#define M_PI 3.14159265358979323846

typedef std::complex<double> complex;

std::vector<complex> Radix2FFT(std::vector<complex> P) {
    //Cooley-Tukey Radix-2 FFT

    int N = P.size();
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
    double theta = (2 * M_PI) / N;
    complex omega_n = std::polar(1.0, theta);
    complex omega = 1;
    for (int i = 0; i < N / 2; i++) {
        res[i] = U_star[i] + omega * V_star[i];
        res[i + N / 2] = U_star[i] - omega * V_star[i];
        omega *= omega_n;
    }
    return res;
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

    // Print the first column values (for verification)
    for (const auto& value : weather_vector) {
        std::cout << value << std::endl;
    }
    return weather_vector;

}


int main() {

    std::string file_path = "../natural_gas_co2_emissions_for_electric_power_sector.csv";
    std::vector<complex> weather_data = get_weather_data(file_path);

    std::vector<complex> result = Radix2FFT(weather_data);
    for (const auto& val : result) {
        std::cout << val << " ";
    }
    std::cout << std::endl;
    return 0;
}