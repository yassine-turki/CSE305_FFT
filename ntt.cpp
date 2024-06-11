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
#include <atomic>
#include <map>
#include "poly_ops.cpp"

////////////////////////////////////////////////////////////////////
// Polynomials with integer coefficents
////////////////////////////////////////////////////////////////////
std::unordered_map<int, std::pair<int, int>> precomputed_roots(){
    std::unordered_map<int, std::pair<int, int>> roots;
     // for (int i = 1; i < 20; i++){
        //     int n = std::pow(2, i);
        //     int p =prime_ntt_find(n);
        //     std::pair<int, int> root = find_2n_roots_parallel(p, n, 10);
        //     roots[n] = root;
        //     std::cout<<i<<" "<<root.first<<" "<<root.second<<std::endl;
        // }
    roots[2] = std::make_pair(2, 3);
    roots[4] = std::make_pair(2, 9);
    roots[8] = std::make_pair(10, 12);
    roots[16] = std::make_pair(19, 46);
    roots[32] = std::make_pair(11, 158);
    roots[64] = std::make_pair(9, 200);
    roots[128] = std::make_pair(3, 86);
    roots[256] = std::make_pair(793, 2838);
    roots[512] = std::make_pair(1254, 49);
    roots[1024] = std::make_pair(11053, 8481);
    roots[2048] = std::make_pair(11054, 11483);
    roots[4096] = std::make_pair(28673, 13657);
    roots[8192] = std::make_pair(81, 8091);
    roots[16384] = std::make_pair(39914, 27768);
    roots[32768] = std::make_pair(3, 21846);
    roots[65536] = std::make_pair(393174, 633878);
    roots[131072] = std::make_pair(357046, 742597);
    return roots;
}

bool is_prime(int p){
    if (p == 0 || p == 1){return false;}
    if (p == 2){return true;}
    int root = static_cast<int>(std::sqrt(p));
    for(int i = 2; i <= root; i++){
        if (p % i == 0){
            return false;
        }
    }
    return true;
}

/* ------------------------------------Parallel is_prime ------------------------------------ */
void is_prime_thread(int p, int start, int end, std::atomic<bool>& prime) {
    for (int i = start; i <= end; ++i) {
        if (p % i == 0) {
            prime = false;
            break;
        }
    }
}

bool is_prime_parallel(int p, size_t num_threads) {
    if (p == 0 || p == 1) return false;
    if (p == 2) return true;

    int root = std::ceil(std::sqrt(p));
    std::atomic<bool> prime(true);

    std::vector<std::thread> threads(num_threads-1);
    int chunk_size = root / num_threads;
    int start_block = 2;
    for (size_t i = 0; i < num_threads - 1; ++i) {
        int end_block = start_block + chunk_size - 1;
        threads[i] = std::thread(is_prime_thread, p, start_block, end_block, std::ref(prime));
        start_block = end_block + 1;
    }
    is_prime_thread(p, start_block, root, prime);

    for (auto& t : threads) {
        t.join();
    }
    return prime;
}


int prime_ntt_find(int n){ 
    // we transformed the algorithm for making sure that the 2n-th root exists
    int p = 2 * n + 1;
    while (true){
        if (is_prime(p)){
            break;
        }
        p += 2 * n;
    }
    return p;
}

std::vector<int> prime_decomposition(int n){
    std::vector<int> prime_factors;
    if (n == 0 or n == 1){return prime_factors;}
    if (is_prime(n)){
        prime_factors.push_back(n);
        return prime_factors;
    }
    for (int i = 2; i <= n / 2 + 1; i++){
        if (is_prime(i) && (n % i == 0)){
            prime_factors.push_back(i);
        }
    }
    return prime_factors;
}

void find_primes_thread(int n, std::vector<int>& prime_factors, int start_block, int end_block, std::mutex& mtx){
    std::vector<int> local_primes;  
    for (int i = start_block; i < end_block; ++i) {
        if (is_prime(i) && (n % i == 0)) {
            local_primes.push_back(i);
        }
    }
    std::lock_guard<std::mutex> lock(mtx);  
    prime_factors.insert(prime_factors.end(), local_primes.begin(), local_primes.end());
}

std::vector<int> prime_decomposition_parallel(int n, int num_threads){
    std::vector<int> prime_factors;
    std::vector<std::thread> threads;
    threads.reserve(num_threads - 1);
    std::mutex mtx;

    if (n == 0 or n == 1){return prime_factors;}
    if (is_prime_parallel(n, num_threads)){
        prime_factors.push_back(n);
        return prime_factors;
    }

    int chunk_size = (n / 2 - 1) / num_threads;

    int start_block = 2;
    for (size_t i = 0; i < num_threads - 1; ++i) {
        int end_block = start_block + chunk_size - 1;
        threads[i] = std::thread(find_primes_thread, n, std::ref(prime_factors), start_block, end_block, std::ref(mtx));
        start_block = end_block;
    }

    find_primes_thread(n, prime_factors, start_block, n / 2 + 1, mtx);

    for (auto& t : threads) {
        t.join();
    }

    return prime_factors;
}

int mod_exp(int base, int exp, int mod){
    int result = 1;
    base = base % mod;
    for (int i = 0; i < exp; i++){
        result = (result * base) % mod;
    }
    if (result < 0){result += mod;}
    return result;
}


/*------------------------------------Find a 2nth root and its inverse------------------------------------*/

std::pair<int, int> find_2n_roots(int p, int n){
    for (int i = 2; i < p; i++){
        int pow = mod_exp(i, n, p);
        if (pow != 1){
            int pow_2 = mod_exp(pow, 2, p);
            if (pow_2 == 1){
                int inv = mod_exp(i, 2 * n - 1, p);
                return {i, inv};
            }
        }
    }
    return {0,0};
}

void find_2n_thread(int p, int n, int start_block, int end_block, std::atomic<int>& generator){
    for (int i = start_block; i <= end_block; i++) {
        int pow = mod_exp(i, n, p);
        if (pow != 1){
            int pow_2 = mod_exp(pow, 2, p);
            if (pow_2 == 1){
                int expected = 0;
                if (generator.compare_exchange_weak(expected, i)){
                    return;
                }
            }
        }
    }
}

std::pair<int, int> find_2n_roots_parallel(int p, int n, size_t num_threads) {
    std::atomic<int> generator(0);
    std::vector<std::thread> threads(num_threads - 1);
    int chunk_size = p / num_threads;
    int start_block = 1;
    for (size_t i = 0; i < num_threads - 1; ++i) {
        int end_block = start_block + chunk_size - 1;
        threads[i] = std::thread(find_2n_thread, p, n, start_block, end_block, std::ref(generator));
        start_block = end_block + 1;
    }

    find_2n_thread(p, n, start_block, p, generator);

    for (auto& t : threads) {
        t.join();
    }

    int root = generator.load();
    int inv_root = mod_exp(root, 2 * n - 1, p); 
    return {root, inv_root};
}
/*------------------------------------Naive Convolution using NTT------------------------------------*/


std::vector<int> ntt(std::vector<int> a, int p, std::unordered_map<int, std::pair<int, int>> roots){
    int n = a.size();
    std::vector<int> a_star(n);
    int psi = roots[n].first;
    for (int i = 0; i < n; i++){
        for (int j = 0; j < n; j++){
            a_star[i] = (a_star[i] + a[j] * mod_exp(psi, 2 * i * j + j, p)) % p;
        }
    }
    return a_star;
}

std::vector<int> intt(std::vector<int> a_star, int p, std::unordered_map<int, std::pair<int, int>> roots){
    int n = a_star.size();
    std::vector<int> a(n);
    int inv_psi = roots[n].second;
    int inv_n = mod_exp(n, p - 2, p);
    for (int i = 0; i < n; i++){
        for(int j = 0; j < n; j++){
            a[i] = (a[i] + a_star[j] * mod_exp(inv_psi, 2 * i * j + i, p)) % p;
        }
        a[i] = (a[i] * inv_n) % p;
    }
    return a;
}

std::vector<int> convolution_ntt(std::vector<int> a, std::vector<int> b, int p, std::unordered_map<int, std::pair<int, int>> roots){
    int n = std::max(a.size(), b.size());
    a.resize(n);
    b.resize(n);
    std::vector<int> a_star = ntt(a, p, roots);
    std::vector<int> b_star = ntt(b, p, roots);
    std::vector<int> c_star(n);
    for (int i = 0; i < n; i++){
        c_star[i] = (a_star[i] * b_star[i]) % p;
    }
    return intt(c_star, p, roots);
}

/*------------------------------------Fast convolution using NTT------------------------------------*/


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

    int psi = find_2n_roots(p, n).first;
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
    std::vector<int> reversed_bit_order_arr(n);
    for (int i = 0; i < n; i++){
        unsigned int pos = reverse_bits(i, num_bits);
        reversed_bit_order_arr[pos] = arr[i];
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

std::vector<int> fast_intt(std::vector<int> a_star, int p, std::unordered_map<int, std::pair<int, int>> roots) {
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

    int psi = roots[n].first;
    int inv_psi = roots[n].second;
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


std::vector<int> convolution_fast(std::vector<int> a, std::vector<int> b, int p, std::unordered_map<int, std::pair<int, int>> roots){
    int n = std::max(a.size(), b.size());
    a.resize(n);
    b.resize(n);
    std::vector<int> a_star = fast_ntt(a, p);
    std::vector<int> b_star = fast_ntt(b, p);
    std::vector<int> c_star(n);
    for (int i = 0; i < n; i++){
        c_star[i] = (a_star[i] * b_star[i]) % p;
    }
    return fast_intt(c_star, p, roots);
}



/*------------------------------------Parallel Naive NTT ------------------------------------*/

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

std::vector<int> ntt_parallel(const std::vector<int>& a, int p, std::unordered_map<int, std::pair<int, int>> roots, size_t num_threads) {
    int n = a.size();
    if (num_threads == 1) {
        return ntt(a, p, roots);
    }

    std::vector<int> a_star(n);
    std::vector<std::thread> threads;
    threads.reserve(num_threads - 1);
    std::mutex mtx;

    int psi = roots[n].first;

    for (size_t i = 1; i < num_threads; ++i) {
        threads.emplace_back(compute_ntt_segment, std::cref(a), std::ref(a_star), p, psi, i, num_threads, n, std::ref(mtx));
    }
    compute_ntt_segment(a, a_star, p, psi, 0, num_threads, n, mtx);

    for (auto& th : threads) {
        th.join();
    }

    return a_star;
}

/*------------------------------------Parallel Naive INTT ------------------------------------*/


void compute_intt_segment(const std::vector<int>& a_star, std::vector<int>& a, size_t start_index, size_t num_threads, int n, int p, int inv_psi, int inv_n, std::mutex& mtx) {
    for (size_t k = start_index; k < n ; k += num_threads) {
        int sum = 0;
        for (size_t j = 0; j < n; ++j) {
            sum = (sum + a_star[j] * mod_exp(inv_psi, 2 * k * j + k, p)) % p;
        }
        std::lock_guard<std::mutex> lock(mtx);
        a[k] = (sum * inv_n) % p;
    }
}

std::vector<int> intt_parallel(std::vector<int>& a_star, int p, std::unordered_map<int, std::pair<int, int>> roots, size_t num_threads) {
    int n = a_star.size();
    if (num_threads == 1) {
        return intt(a_star, p, roots);
    }

    std::vector<std::thread> threads;
    threads.reserve(num_threads - 1);

    std::vector<int> a(n);
    std::mutex mtx;

    int inv_psi = roots[n].second;
    int inv_n = mod_exp(n, p - 2, p);

    for (size_t i = 1; i < num_threads; ++i) {
        threads.emplace_back(compute_intt_segment, std::cref(a_star), std::ref(a), i, num_threads, n, p, inv_psi, inv_n, std::ref(mtx));
    }
    compute_intt_segment(a_star, a, 0, num_threads, n, p, inv_psi, inv_n, mtx);

    for (auto& th : threads) {
        th.join();
    }

    return a;
}

/*------------------------------------Parallel naive convolution ------------------------------------*/

void compute_convolution_segment(const std::vector<int>& a, const std::vector<int>& b, std::vector<int>& c, size_t start_index, size_t num_threads, int n, int p, std::mutex& mtx){
    for (size_t k = start_index; k < n ; k += num_threads) {
        c[k] = (a[k] * b[k]) % p;
    }
}

std::vector<int> convolution_ntt_parallel(std::vector<int> a, std::vector<int> b, int p, std::unordered_map<int, std::pair<int, int>> roots, int num_threads){
    int n = std::max(a.size(), b.size());
    a.resize(n);
    b.resize(n);
    std::vector<int> a_star = ntt_parallel(a, p, roots, num_threads / 2);
    std::vector<int> b_star = ntt_parallel(b, p, roots, num_threads - num_threads / 2);
    std::vector<int> c_star(n);

    std::vector<std::thread> threads;
    threads.reserve(num_threads - 1);
    std::mutex mtx;

    for (size_t i = 1; i < num_threads; ++i) {
        threads.emplace_back(compute_convolution_segment, std::cref(a_star), std::cref(b_star), std::ref(c_star), i, num_threads, n, p, std::ref(mtx));
    }
    compute_convolution_segment(a_star, b_star, c_star, 0, num_threads, n, p, mtx);

    for (auto& th : threads) {
        th.join();
    }

    return intt_parallel(c_star, p, roots, num_threads);
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

    int psi = find_2n_roots(p, n).first;
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
    std::thread t2(fast_ntt_parallel, std::ref(v), p, num_threads - num_threads / 2);
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


/*------------------------------------Parallel INTT Gentleman Sande ------------------------------------*/

void compute_reverse_bit_order_segment(const std::vector<int>& arr, std::vector<int>& reversed_bit_order_arr, size_t start_index, size_t num_threads, int n, int num_bits, std::mutex& mtx){
    for (size_t k = start_index; k < n ; k += num_threads) {
        unsigned int pos = reverse_bits(k, num_bits);
        std::lock_guard<std::mutex> lock(mtx);
        reversed_bit_order_arr[pos] = arr[k];
    }
}

std::vector<int> reverse_bit_order_array_parallel(std::vector<int>& arr, int num_threads){
    int n = arr.size();
    int num_bits = std::log2(n);
    std::vector<int> reversed_bit_order_arr(n);
    std::vector<std::thread> threads;
    threads.reserve(num_threads - 1);
    std::mutex mtx;

    for (size_t i = 1; i < num_threads; ++i) {
        threads.emplace_back(compute_reverse_bit_order_segment, std::cref(arr), std::ref(reversed_bit_order_arr), i, num_threads, n, num_bits, std::ref(mtx));
    }
    compute_reverse_bit_order_segment(arr, reversed_bit_order_arr, 0, num_threads, n, num_bits, mtx);

    for (auto& th : threads) {
        th.join();
    }

    return reversed_bit_order_arr;
}

std::vector<int> compute_powers_of_psi_parallel(int psi, int n, int p, int num_threads) {
    int num_bits = std::log2(n);
    std::vector<int> powers(n);
    std::vector<int> bit_reversed_powers(n);

    powers[0] = 1;
    for (int i = 1; i < n; i++) {
        powers[i] = (powers[i - 1] * psi) % p;
    }

    std::vector<std::thread> threads;
    threads.reserve(num_threads - 1);
    std::mutex mtx;

    for (size_t i = 1; i < num_threads; ++i) {
        threads.emplace_back(compute_reverse_bit_order_segment, std::cref(powers), std::ref(bit_reversed_powers), i, num_threads, n, num_bits, std::ref(mtx));
    }
    compute_reverse_bit_order_segment(powers, bit_reversed_powers, 0, num_threads, n, num_bits, mtx);

    for (auto& th : threads) {
        th.join();
    }

    return bit_reversed_powers;
}

void compute_segment_fast_intt(int start, int end, int t, int p, const std::vector<int>& powers_inv_psi, std::vector<int>& a_star, std::mutex& mtx) {
    for (int i = start; i < end; i++) {
        int S = powers_inv_psi[end + i];
        int k = 2 * t * i;
        for (int j = k; j < k + t; j++) {
            std::lock_guard<std::mutex> lock(mtx);
            int U = a_star[j];
            int V = a_star[j + t];
            a_star[j] = (U + V) % p;
            if (a_star[j] < 0) a_star[j] += p;
            int W = (U - V) % p;
            if (W < 0) W += p;
            a_star[j + t] = (W * S) % p;
            if (a_star[j + t] < 0) a_star[j + t] += p;
        }
    }
}

std::vector<int> fast_intt_parallel(std::vector<int> a_star, int p, std::unordered_map<int, std::pair<int, int>> roots, int num_threads) {
    int n = a_star.size();
    if (n == 1) {
        return a_star;
    }

    if (!is_power_of_two(n)) {
        std::cout << "The input size " << n << " is not a power of 2! Resizing input" << std::endl;
        a_star.resize(next_power_of_two(n));
        n = a_star.size();
    }

    a_star = reverse_bit_order_array_parallel(a_star, num_threads);

    int psi = roots[n].first;
    int inv_psi = roots[n].second;
    int inv_n = mod_exp(n, p - 2, p);
    std::vector<int> powers_inv_psi = compute_powers_of_psi_parallel(inv_psi, n, p, num_threads);
    
    int t = 1;
    int m = n / 2;
    std::mutex mtx;

    while (m > 0) {
        std::vector<std::thread> threads;
        int segment_size = m / num_threads;

        for (int i = 0; i < num_threads; i++) {
            int start = i * segment_size;
            int end = (i == num_threads - 1) ? m : start + segment_size;
            threads.emplace_back(compute_segment_fast_intt, start, end, t, p, std::ref(powers_inv_psi), std::ref(a_star), std::ref(mtx));
        }

        for (auto& thread : threads) {
            thread.join();
        }

        t *= 2;
        m /= 2;
    }

    std::vector<std::thread> scaling_threads;
    int segment_size = n / num_threads;

    for (int i = 0; i < num_threads; i++) {
        int start = i * segment_size;
        int end;
        if (i == num_threads - 1){
            end = n;
        }
        else {
            end = start + segment_size;
        } 
        scaling_threads.emplace_back([start, end, inv_n, p, &a_star]() {
            for (int i = start; i < end; i++) {
                a_star[i] = (a_star[i] * inv_n) % p;
                if (a_star[i] < 0) a_star[i] += p;
            }
        });
    }

    for (auto& thread : scaling_threads) {
        thread.join();
    }

    return a_star;
}


/*------------------------------------Parallel fast convolution ------------------------------------*/

std::vector<int> convolution_fast_parallel(std::vector<int> a, std::vector<int> b, int p, std::unordered_map<int, std::pair<int, int>> roots, int num_threads){
    int n = std::max(a.size(), b.size());
    a.resize(n);
    b.resize(n);
    fast_ntt_parallel(a, p, num_threads / 2);
    fast_ntt_parallel(b, p, num_threads - num_threads / 2);
    std::vector<int> c(n);

    std::vector<std::thread> threads;
    threads.reserve(num_threads - 1);
    std::mutex mtx;

    for (size_t i = 1; i < num_threads; ++i) {
        threads.emplace_back(compute_convolution_segment, std::cref(a), std::cref(b), std::ref(c), i, num_threads, n, p, std::ref(mtx));
    }
    compute_convolution_segment(a, b, c, 0, num_threads, n, p, mtx);

    for (auto& th : threads) {
        th.join();
    }
    return fast_intt_parallel(c, p, roots, num_threads);
}