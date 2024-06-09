#include "plot.cpp"

void compress_block_threads(
        std::vector<complex> &block, double cut_off, int start, int end
) {
    for (int i = start; i < end; i++) {
        if (abs(block[i]) < cut_off){
            block[i] = 0;
        }
    }
}

void compress(std::vector<complex> &data, double cut_off, size_t num_threads) {
    int N = data.size();
    int block_size = N/num_threads;
    std::vector<std::thread> threads(num_threads - 1);
    int start = 0;
    for (int i = 0; i < (num_threads - 1); i++) {
        int end = start + block_size;
        threads[i] = std::thread(compress_block_threads,std::ref(data), cut_off, start, end);
        start += block_size;
    }
    compress_block_threads(data, cut_off, start, N);
    for(auto &th : threads) {
        th.join();
    }
}