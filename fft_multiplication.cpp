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
    for (int i = 1; i < p; i++){
        int ct = 0;
        for (int prime: prime_factors){
            int y = mod_exp(i, (p - 1) / prime, p);
            if (y == 1){
                break;
            }
            ct += 1;
            if (ct == prime_factors.size()){
                return i;
            }
        }
    }
    return 0;
}


int find_omega(int p, int n){
    int g = find_generator(p, n);
    int k = (p - 1) / n;
    return mod_exp(g, k, p);
}

std::vector<int> find_2n_roots(int p, int n){
    int g = find_generator(p, n);
    int k = (p - 1) / n;
    if (k % 2 != 0){
        throw std::runtime_error("There is no 2n-th root of unity!");
    }
    if (g == 0){
        throw std::runtime_error("There is no generator!");
    }
    int root_1 = mod_exp(g, k / 2, p);
    int root_2 = mod_exp(root_1, 2 * n - 1, p);
    std::vector<int> roots = {root_1, root_2};
    return roots;
}

std::vector<int> ntt(std::vector<int> a, int p){
    int n = a.size();
    std::vector<int> a_star(n);
    int psi = find_2n_roots(p, n)[0];
    std::cout<<"Psi in ntt is "<<psi<<std::endl;
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
    int inv_psi = find_2n_roots(p, n)[1];
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

    int psi = find_2n_roots(p, n)[0];
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

unsigned int reverse_bits(unsigned int x, int num_bits) {
    unsigned int result = 0;
    for (int i = 0; i < num_bits; i++) {
        result = (result << 1) | (x & 1);
        x >>= 1;
    }
    return result;
}

std::vector<int> reverse_bit_order_array(std::vector<int> arr){
    int n = arr.size();
    int num_bits = std::log2(n);
    std::vector<unsigned int> new_positions(n);
    for (unsigned int i = 0; i < n; i++){
        unsigned int pos = reverse_bits(i, num_bits);
        new_positions[i] = pos;
    }
    std::vector<int> reversed_bit_order_arr(n);
    for (int i = 0; i < n; i++){
        unsigned int new_pos = new_positions[i];
        reversed_bit_order_arr[new_pos] = arr[i];
    }
    return reversed_bit_order_arr;
}

std::vector<int> compute_powers_of_psi(int psi, int n, int p) {
    int num_bits = std::log2(n);
    std::vector<int> powers(n);
    std::vector<int> bit_reversed_powers(n);

    powers[0] = 1;
    for (int i = 1; i < n; i++) {
        powers[i] = (powers[i - 1] * psi) % p;
    }

    for (int i = 0; i < n; i++) {
        unsigned int reversed_index = reverse_bits(i, num_bits);
        bit_reversed_powers[reversed_index] = powers[i];
    }

    return bit_reversed_powers;
}

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

    a_star = reverse_bit_order_array(a_star);

    std::vector<int> roots = find_2n_roots(p, n);
    int psi = roots[0];
    int inv_psi = roots[1];
    if (((inv_psi * psi) % p) != 1){
        throw std::runtime_error("The roots are not well computed");
    }
    int inv_n = mod_exp(n, p - 2, p);
    std::vector<int> powers_inv_psi = compute_powers_of_psi(inv_psi, n, p);
    
    int t = 1;
    int m = n / 2;
    while (m > 0) {
        int k = 0;
        for (int i = 0; i < m; i++) {
            int S = powers_inv_psi[m + i];
            for (int j = k; j < k + t; j++) {
                int U = a_star[j];
                int V = a_star[j + t];
                a_star[j] = (U + V) % p;
                if (a_star[j] < 0) a_star[j] += p;
                int W = (U - V) % p;
                if (W < 0) W += p;
                a_star[j + t] = (W * S) % p;
                if (a_star[j + t] < 0) a_star[j + t] += p;
            }
            k += 2 * t;
        }
        t *= 2;
        m /= 2;
    }
    for (int i = 0; i < n; i++) {
        a_star[i] = (a_star[i] * inv_n) % p;
        if (a_star[i] < 0) a_star[i] += p;
    }
    return a_star;
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



/*------------------------------------Parallel NTT ------------------------------------*/

void compute_ntt_segment(const std::vector<int>& a, std::vector<int>& a_star, int p, int psi, size_t start_index, size_t num_threads, int n, std::mutex& mtx) {
    for (size_t i = start_index; i < n; i += num_threads) {
        int sum = 0;
        for (size_t j = 0; j < n; ++j) {
            sum = (sum + a[j] * mod_exp(psi, 2 * i * j + j, p)) % p;
        }
        std::lock_guard<std::mutex> lock(mtx);
        a_star[i] = sum;
    }
}

std::vector<int> ntt_parallel(const std::vector<int>& a, int p, size_t num_threads) {
    int n = a.size();
    if (num_threads == 1) {
        return ntt(a, p);
    }

    std::vector<int> a_star(n);
    std::vector<std::thread> threads;
    threads.reserve(num_threads - 1);
    std::mutex mtx;

    int psi = find_2n_roots(p, n)[0];

    for (size_t i = 1; i < num_threads; ++i) {
        threads.emplace_back(compute_ntt_segment, std::cref(a), std::ref(a_star), p, psi, i, num_threads, n, std::ref(mtx));
    }
    compute_ntt_segment(a, a_star, p, psi, 0, num_threads, n, mtx);

    for (auto& th : threads) {
        th.join();
    }

    return a_star;
}

/*------------------------------------Parallel NTT Cooley ------------------------------------*/

void split_data_NTT(const std::vector<int>& a, std::vector<int>& u, std::vector<int>& v, size_t start, size_t end) {
    for (size_t i = start; i < end; ++i) {
        u[i] = a[2 * i];
        v[i] = a[2 * i + 1];
    }
}

void combine_results_NTT(std::vector<int>& a, const std::vector<int>& u, const std::vector<int>& v, int p, int psi, size_t n, size_t start, size_t end) {
    int t = mod_exp(psi, 2 * start+1, p);
    int psi_power = mod_exp(psi, 2, p); // psi squared
    for (size_t i = start; i < end; ++i) {
        a[i] = (u[i] + t * v[i]) % p;
        if (a[i] < 0) a[i] += p;
        a[i + n / 2] = (u[i] - t * v[i]) % p;
        if (a[i + n / 2] < 0) a[i + n / 2] += p;
        t = (t * psi_power) % p;
    }
}

void fast_ntt_parallel(std::vector<int>& a, int p, size_t num_threads) {
    size_t n = a.size();
    if (!is_power_of_two(n)) {
        std::cout << "The input size " << n << " is not a power of 2! Resizing input" << std::endl;
        n = next_power_of_two(n);
        a.resize(n);
    }

    if (num_threads == 1 || n == 1) {
        a = fast_ntt(a, p);
        return;
    }

    int psi = find_2n_roots(p, n)[0];
    std::vector<int> u(n / 2), v(n / 2);

    size_t split_threads = std::min(num_threads, n / 2);
    size_t chunk_size = (n / 2 + split_threads - 1) / split_threads;
    std::vector<std::thread> split_threads_list;

    for (size_t t = 0; t < split_threads; ++t) {
        size_t start = t * chunk_size;
        size_t end = std::min(start + chunk_size, n / 2);
        split_threads_list.emplace_back(split_data_NTT, std::cref(a), std::ref(u), std::ref(v), start, end);
    }

    for (auto& th : split_threads_list) {
        th.join();
    }
    std::thread t1(fast_ntt_parallel, std::ref(u), p, num_threads / 2);
    std::thread t2(fast_ntt_parallel, std::ref(v), p, num_threads / 2);
    t1.join();
    t2.join();

    split_threads_list.clear();
    for (size_t t = 0; t < split_threads; ++t) {
        size_t start = t * chunk_size;
        size_t end = std::min(start + chunk_size, n / 2);
        split_threads_list.emplace_back(combine_results_NTT, std::ref(a), std::cref(u), std::cref(v), p, psi, n, start, end);
    }

    for (auto& th : split_threads_list) {
        th.join();
    }
}
