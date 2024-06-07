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
#include <set>
#include <random>
#include "fft.cpp"


typedef std::complex<double> complex;

template<typename T>
T evaluate_polynomial_trivial(std::vector<T> P, T x){
    int N = P.size();
    T X = 1;
    T y;
    for (int j = 0; j < N; j++){
        y += P[j] * X;
        X *= x;
    }
    return y;
}

template<typename T>
std::vector<T> add_polynomials(std::vector<T> P, std::vector<T> Q){
    int N = std::max(P.size(), Q.size());
    P.resize(N);
    Q.resize(N);
    std::vector<T> R(N);
    for (int j = 0; j < N; j++){
        R[j] = P[j] + Q[j];
    }
    return R;
}

template<typename T>
std::vector<T> multiply_polynomials_trivial(std::vector<T> P, std::vector<T> Q){
    int N = P.size();
    int M = Q.size();
    std::vector<T> R(M + N, 0.);
    for (int j = 0; j < N; j++){
        for (int k = 0; k < M; k++){
            R[j + k] += P[j] * Q[k];
        }
    }
    return R;
}

template<typename T>
T horner_evaluate(std::vector<T> P, T x){
    int N = P.size();
    T y = P[N - 1];
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

////////////////////////////////////////////////////////////////////
// Polynomials with integer coefficents
////////////////////////////////////////////////////////////////////

bool is_prime(int p){
    if (p == 0 || p == 1){return false;}
    if (p == 2){return true;}
    int root = std::ceil(std::sqrt(p));
    for(int i = 2; i <= root; i++){
        if (p % i == 0){
            return false;
        }
    }
    return true;
}

std::vector<int> prime_fft_find(int n){ 
    // we transformed the algorithm for making sure that the 2n-th root exists
    int k = 2;
    int p = n + 1;
    while (true){
        if (is_prime(p)){
            break;
        }
        k += 2;
        p += 2 * n;
    }
    std::vector<int> prime_pair = {p, k};
    return prime_pair;
}

std::vector<int> prime_decomposition(int n){
    std::vector<int> prime_factors;
    if (n == 0 or n == 1){return prime_factors;}
    if (is_prime(n)){
        prime_factors.push_back(n);
        return prime_factors;
    }
    for (int i = 2; i <= n; i++){
        if (is_prime(i) && (n % i == 0)){
            prime_factors.push_back(i);
        }
    }
    return prime_factors;
}

int mod_exp(int base, int exp, int mod) {
    int result = 1;
    base = base % mod;
    while (exp > 0) {
        if (exp % 2 == 1) {
            result = (result * base) % mod;
        }
        exp = exp >> 1; 
        base = (base * base) % mod; 
    }

    return result;
}

int find_generator(int p, int n){
    int k = (p - 1) / n;
    std::vector<int> prime_factors = prime_decomposition(k);
    auto it = std::find(prime_factors.begin(), prime_factors.end(), 2);
    if (it == prime_factors.end()){prime_factors.push_back(2);}
    bool found = false;
    int g = 0;
    while(!found){
        g = rand() % (p - 1) + 1;
        found = true;
        for (int prime: prime_factors){
            int y = mod_exp(g, (p - 1) / prime, p);
            if (y == 1){
                found = false;
            }
        }
    }
    return g;
}

int find_omega(int p, int n){
    int g = find_generator(p, n);
    int k = (p - 1) / n;
    return mod_exp(g, k, p);
}

int find_2n_root(int p, int n){
    int g = find_generator(p, n);
    int k = (p - 1) / n;
    if (k % 2 != 0){
        throw std::runtime_error("There is no 2n-th root of unity!");
    }
    return mod_exp(g, k / 2, p);
}

std::vector<int> ntt(std::vector<int> a, int p){
    int n = a.size();
    std::vector<int> a_star(n);
    // std::vector<int> pair = prime_fft_find(n);
    // int p = pair[0];
    int psi = find_2n_root(p, n);
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            a_star[i] = (a_star[i] + a[j] * mod_exp(psi, 2 * i * j + j, p)) % p;
        }
    }
    return a_star;
}

std::vector<int> intt(std::vector<int> a_star, int p){
    int n = a_star.size();
    std::vector<int> a(n);
    // std::vector<int> pair = prime_fft_find(n);
    // int p = pair[0];
    int psi = find_2n_root(p, n);
    int inv_psi = mod_exp(psi, 2 * n - 1, p);
    int inv_n = mod_exp(n, p - 2, p);
    for (int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            a[i] = (a[i] + a_star[j] * mod_exp(inv_psi, 2 * i * j + i, p)) % p;
        }
        a[i] = (a[i] * inv_n) % p;
    }
    return a;
}

std::vector<int> convolution_ntt(std::vector<int> a, std::vector<int> b){
    int n = std::max(a.size(), b.size());
    a.resize(n);
    b.resize(n);
    std::vector<int> pair = prime_fft_find(n);
    int p = pair[0];
    std::vector<int> a_star = ntt(a, p);
    std::vector<int> b_star = ntt(b, p);
    std::vector<int> c_star(n);
    for (int i = 0; i < n; i++){
        c_star[i] = (a_star[i] * b_star[i]) % p;
    }
    return intt(c_star, p);
}

std::vector<int> Radix2FFT_int(std::vector<int> P, int p, int omega) {
    int N = P.size();
    if(!is_power_of_two(N)){
        std::cout << "The input size " << N << " is not a power of 2! Resizing input" << std::endl;
        P.resize(next_power_of_two(P.size()));
        N = P.size();
    }

    if (N == 1) {
        return P;
    }

    // std::vector<int> pair = prime_fft_find(N);
    // int p = pair[0];
    // int k = pair[1];
    // int g = find_generator(p, N);
    // int omega = mod_exp(g, k, p);

    std::vector<int> U(N / 2), V(N / 2);
    for (int i = 0; i < N / 2; i++) {
        U[i] = P[2 * i];
        V[i] = P[2 * i + 1];
    }

    std::vector<int> U_star = Radix2FFT_int(U, p, mod_exp(omega, 2, p));
    std::vector<int> V_star = Radix2FFT_int(V, p, mod_exp(omega, 2, p));
    std::vector<int> res(N);

    int t = 1;
    for (int i = 0; i < N / 2; i++) {
        res[i] = U_star[i] + t * V_star[i];
        res[i + N / 2] = U_star[i] - t * V_star[i];
        t = (t * omega) % p;
    }
    return res;
}

std::vector<int> InverseRadix2FFT_int(std::vector<int> P_star, int p, int omega) {
    int N = P_star.size();
    if (N == 1) {
        return P_star;
    }

    if (!is_power_of_two(N)) {
        std::cout << "The input size " << N << " is not a power of 2! Resizing input" << std::endl;
        P_star.resize(next_power_of_two(N));
        N = P_star.size();
    }

    // std::vector<int> pair = prime_fft_find(N);
    // int p = pair[0];
    // int k = pair[1];
    // int g = find_generator(p, N);
    // int omega = mod_exp(g, k, p);

    std::vector<int> U_star(N / 2), V_star(N / 2);

    for (int j = 0; j < N / 2; j++) {
        U_star[j] = P_star[2 * j];
        V_star[j] = P_star[2 * j + 1];
    }

    std::vector<int> U = InverseRadix2FFT_int(U_star, p, mod_exp(omega, 2, p));
    std::vector<int> V = InverseRadix2FFT_int(V_star, p, mod_exp(omega, 2, p));

    int t = 1;
    std::vector<int> P(N);

    for (int j = 0; j < N / 2; j++) {
        P[j] = (U[j] + t * V[j]) / 2;
        P[j + N / 2] = (U[j] - t * V[j]) / 2;
        t = (t * omega) % p;
    }

    return P;
}

std::vector<int> multiply_polynomials_FFT_int(std::vector<int> P, std::vector<int> Q){
    int N = P.size();
    int M = Q.size();
    int l = next_power_of_two(M + N);
    P.resize(l);
    Q.resize(l);

    std::vector<int> pair = prime_fft_find(l);
    int p = pair[0];  
    int k = pair[1];
    int g = find_generator(p, N);
    int omega = mod_exp(g, k, p);

    std::vector<int> P_star = Radix2FFT_int(P, p, omega);
    std::vector<int> Q_star = Radix2FFT_int(Q, p, omega);
    std::vector<int> R_star(l);
    for (int j = 0; j < l; j++){
        R_star[j] = (P_star[j] * Q_star[j]) % p;
    }
    return InverseRadix2FFT_int(R_star, p, mod_exp(omega, l - 1, p));
}

