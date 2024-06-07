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
    int p = 2 * n + 1;
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

std::vector<int> fast_ntt(std::vector<int> a, int p) {
    int n = a.size();
    if(!is_power_of_two(n)){
        std::cout << "The input size " << n << " is not a power of 2! Resizing input" << std::endl;
        a.resize(next_power_of_two(a.size()));
        n = a.size();
    }

    if (n == 1) {
        return a;
    }

    int psi = find_2n_root(p, n);
    std::vector<int> u(n / 2), v(n / 2);
    for (int i = 0; i < n / 2; i++) {
        u[i] = a[2 * i];
        v[i] = a[2 * i + 1];
    }

    std::vector<int> u_star = fast_ntt(u, p);
    std::vector<int> v_star = fast_ntt(v, p);
    std::vector<int> a_star(n);

    int t = psi;
    for (int i = 0; i < n / 2; i++) {
        a_star[i] = (u_star[i] + t * v_star[i]) % p;
        if (a_star[i] < 0) a_star[i] += p;
        a_star[i + n / 2] = (u_star[i] - t * v_star[i]) % p;
        if (a_star[i + n / 2] < 0) a_star[i + n / 2] += p;
        t = (t * mod_exp(psi, 2, p)) % p;
    }
    return a_star;
}

// std::vector<int> fast_ntt_2(std::vector<int> a, int p){
//     int n = a.size();
//     if(!is_power_of_two(n)){
//         std::cout << "The input size " << n << " is not a power of 2! Resizing input" << std::endl;
//         a.resize(next_power_of_two(a.size()));
//         n = a.size();
//     }

//     if (n == 1) {
//         return a;
//     }

//     int psi = find_2n_root(p, n);

//     int t = n;
//     for (int m = 1; m < n; m = 2 * m){
//         t = t / 2;
//         for (int i = 0; i < m; i++){
//             int j_1 = 2 * i * t;
//             int j_2 = j_1 + t - 1;
//             int S = mod_exp(psi, 2 * n - m - i, p);
//             for (int j = j_1; j < j_2; j++){
//                 int u = a[j];
//                 int v = a[j + t] * S;
//                 a[j] = (u + v) % p;
//                 a[j + t] = (u - v) % p;
//             }
//         }
//     }
//     return a;
// }


std::vector<int> fast_intt(std::vector<int> a_star, int p) {
    int n = a_star.size();
    if (n == 1) {
        return a_star;
    }

    if (!is_power_of_two(n)) {
        std::cout << "The input size " << n << " is not a power of 2! Resizing input" << std::endl;
        a_star.resize(next_power_of_two(n));
        n = a_star.size();
    }

    int psi = find_2n_root(p, n);
    int inv_psi = mod_exp(psi, 2 * n - 1, p);
    int inv_n = mod_exp(n, p - 2, p);

    std::vector<int> u_star(n / 2), v_star(n / 2);

    for (int j = 0; j < n / 2; j++) {
        u_star[j] = a_star[j];
        v_star[j] = a_star[n / 2 + j];
    }

    std::vector<int> u = fast_intt(u_star, p);
    std::vector<int> v = fast_intt(v_star, p);

    int t = 1;
    std::vector<int> a(n);

    for (int i = 0; i < n / 2; i++) {
        a[2 * i] = ((u[i] + v[i]) * inv_psi) % p;
        a[2 * i] = (a[2 * i] * inv_n) % p;
        a[2 * i + 1] = ((u[i] - v[i]) * inv_psi) % p;
        a[2 * i + 1] = (a[2 * i + 1] * inv_n) % p;
    }

    return a;
}

std::vector<int> FFT_convolution(std::vector<int> a, std::vector<int> b){
    int n = std::max(a.size(), b.size());
    a.resize(n);
    b.resize(n);
    std::vector<int> pair = prime_fft_find(n);
    int p = pair[0];
    std::vector<int> a_star = fast_ntt(a, p);
    std::vector<int> b_star = fast_ntt(b, p);
    std::vector<int> c_star(n);
    for (int i = 0; i < n; i++){
        c_star[i] = (a_star[i] * b_star[i]) % p;
    }
    return fast_intt(c_star, p);
}
