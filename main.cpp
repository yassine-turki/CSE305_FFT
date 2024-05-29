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
    if (N == 1){
        return P_star;
    }

    std::vector<complex> U_star(N / 2), V_star(N / 2);

    for (int j = 0; j < N / 2; j++){
        U_star[j] = P_star[2 * j];
        V_star[j] = P_star[2 * j + 1];
    }

    std::vector<complex> U = InverseRadix2FFT(U_star);
    std::vector<complex> V = InverseRadix2FFT(V_star);

    double theta = 2 * M_PI - (2 * M_PI) / N;
    complex omega_n = std::polar(1.0, theta);
    complex omega = 1;
    std::vector<complex> P(N);

    for (int j = 0; j < N / 2; j++){
        P[j] = (U[j] + omega * V[j]) / 2.;
        P[j + N / 2] = (U[j] - omega * V[j]) / 2.;
        omega = omega * omega_n;
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
    return weather_vector;

}


int main() {

//    std::string file_path = "../data/natural_gas_co2_emissions_for_electric_power_sector.csv";
//    std::string file_path = "../data/Paris_data.csv";
//    std::vector<complex> weather_data = get_weather_data(file_path);
std::vector<complex> weather_data = {1, 2, 3, 4};

    std::vector<complex> result = Radix2FFT(weather_data);
    for (int i = 0; i < result.size(); i++) {
        std::cout << result[i] << std::endl;
    }

//    std::vector<complex> inverse_result = InverseRadix2FFT(result);
//    if(weather_data.size() != inverse_result.size()) {
//        std::cout << "Error: The size of the input and output vectors are not the same!" << std::endl;
//    }
//
//    for(int i = 0; i < inverse_result.size(); i++) {
//        if(weather_data[i] != inverse_result[i]) {
//            std::cout << "Error at index " << i << std::endl;
//            std::cout << "Input: " << weather_data[i] << " Output: " << inverse_result[i] << std::endl;
//        }

    return 0;
}