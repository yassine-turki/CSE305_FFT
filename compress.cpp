#include "data.cpp"
#define STB_IMAGE_IMPLEMENTATION
#include "STB/stb_image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "STB/stb_image_write.h"



/*------------------------ Image FFT ------------------------*/


// Function to read an image and convert to grayscale
std::vector<complex>* read_img(const std::string& filename, int &width, int &height) {
    int channels;
    uint8_t* rgb_img = stbi_load(filename.c_str(), &width, &height, &channels, 3);
    if (rgb_img == nullptr) {
        std::cerr << "Error loading image: " << filename << std::endl;
        return nullptr;
    }

    std::cout << "Image Size: " << width << "x" << height << " (" << channels << " channels)" << std::endl;

    size_t grey_img_size = width * height;
    std::vector<complex> img_data(grey_img_size);

    for (int i = 0; i < grey_img_size; ++i) { // Convert to grey scale (formula from https://en.wikipedia.org/wiki/Grayscale)
        uint8_t* pixel = rgb_img + 3 * i;
        double grey_value = 0.299 * pixel[0] + 0.587 * pixel[1] + 0.114 * pixel[2];
        img_data[i] = grey_value;
    }

    stbi_image_free(rgb_img);

    std::cout << "Image read successfully!" << std::endl;
    return new std::vector<complex>(std::move(img_data));
}

// Function to write an image
void write_img(const std::vector<complex> &img_data, int width, int height, const std::string& filename) {
    std::cout << "Writing image to: " << filename << std::endl;

    std::vector<uint8_t> out_img(width * height);

    for (size_t i = 0; i < img_data.size(); ++i) {
        double value = std::max(0.0, std::min(255.0, img_data[i].real()));
        out_img[i] = static_cast<uint8_t>(value);
    }

    stbi_write_png(filename.c_str(), width, height, 1, out_img.data(), width);
}



// Function to compute FFT for rows
void compute_fft_row(std::vector<complex>& img, int width, int height, int start, int end, int num_sub_threads = 1) {
    for (int i = start; i < end; ++i) {
        std::vector<complex> row(img.begin() + i * width, img.begin() + (i + 1) * width);
        if (num_sub_threads == 1) {
            Radix2FFT(row);
        } else {
            FFT_Radix2_Parallel(row, num_sub_threads);
        }
        std::copy(row.begin(), row.end(), img.begin() + i * width);
    }
}

// Function to compute FFT for columns
void compute_fft_column(std::vector<complex>& img, int width, int height, int start, int end, int num_sub_threads = 1) {
    for (int j = start; j < end; ++j) {
        std::vector<complex> column(height);
        std::copy(&img[j], &img[j + height * width], column.begin());

        if (num_sub_threads == 1) {
            Radix2FFT(column);
        } else {
            FFT_Radix2_Parallel(column, num_sub_threads);
        }

        std::copy(column.begin(), column.end(), &img[j]);
    }
}
//void compute_fft_column(std::vector<complex>& img, int width, int height, int start, int end, int num_sub_threads = 1) {
//    for (int j = start; j < end; ++j) {
//        std::vector<complex> column(height);
//        for (int i = 0; i < height; ++i) {
//            column[i] = img[i * width + j];
//        }
//
//        if (num_sub_threads == 1) {
//            Radix2FFT(column);
//        } else {
//            FFT_Radix2_Parallel(column, num_sub_threads);
//        }
//
//        for (int i = 0; i < height; ++i) {
//            img[i * width + j] = column[i];
//        }
//    }
//}

// Function to compute inverse FFT for rows
void compute_ifft_row(std::vector<complex>& img, int width, int height, int start, int end, int num_sub_threads = 1) {
    for (int i = start; i < end; ++i) {
        std::vector<complex> row(img.begin() + i * width, img.begin() + (i + 1) * width);

        if (num_sub_threads == 1) {
            InverseRadix2FFT(row);
        } else {
            Inverse_Radix2_Parallel(row, num_sub_threads);
        }

        std::copy(row.begin(), row.end(), img.begin() + i * width);
    }
}

// Function to compute inverse FFT for columns
void compute_ifft_column(std::vector<complex>& img, int width, int height, int start, int end, int num_sub_threads = 1) {
    for (int j = start; j < end; ++j) {
        std::vector<complex> column(height);
        for (int i = 0; i < height; ++i) {
            column[i] = img[i * width + j];
        }

        if (num_sub_threads == 1) {
            InverseRadix2FFT(column);
        } else {
            Inverse_Radix2_Parallel(column, num_sub_threads);
        }

        for (int i = 0; i < height; ++i) {
            img[i * width + j] = column[i];
        }
    }
}

// Function to compute FFT of an image in parallel
void parallel_fft_image(std::vector<complex>& img, int width, int height, int num_threads, bool double_parallel = false) {
    int num_sub_threads = 1;
    if (double_parallel) {
        num_threads /= 2;
        num_sub_threads = 2;
    }

    int block_size = height / num_threads;
    std::vector<std::thread> threads;

    for (int i = 0; i < num_threads - 1; ++i) {
        int start = i * block_size;
        int end = (i + 1) * block_size;
        threads.emplace_back(compute_fft_row, std::ref(img), width, height, start, end, num_sub_threads);
    }
    compute_fft_row(img, width, height, (num_threads - 1) * block_size, height, num_sub_threads);

    for (auto& th : threads) {
        th.join();
    }
    threads.clear();

    block_size = width / num_threads;
    for (int i = 0; i < num_threads - 1; ++i) {
        int start = i * block_size;
        int end = (i + 1) * block_size;
        threads.emplace_back(compute_fft_column, std::ref(img), width, height, start, end, num_sub_threads);
    }
    compute_fft_column(img, width, height, (num_threads - 1) * block_size, width, num_sub_threads);

    for (auto& th : threads) {
        th.join();
    }
}

// Function to compute inverse FFT of an image in parallel
void parallel_inverse_fft_image(std::vector<complex>& img, int width, int height, int num_threads, bool double_parallel = false) {
    int num_sub_threads = 1;
    if (double_parallel) {
        num_threads /= 2;
        num_sub_threads = 2;
    }

    int block_size = width / num_threads;
    std::vector<std::thread> threads;

    for (int i = 0; i < num_threads - 1; ++i) {
        int start = i * block_size;
        int end = (i + 1) * block_size;
        threads.emplace_back(compute_ifft_column, std::ref(img), width, height, start, end, num_sub_threads);
    }
    compute_ifft_column(img, width, height, (num_threads - 1) * block_size, width, num_sub_threads);

    for (auto& th : threads) {
        th.join();
    }
    threads.clear();

    block_size = height / num_threads;
    for (int i = 0; i < num_threads - 1; ++i) {
        int start = i * block_size;
        int end = (i + 1) * block_size;
        threads.emplace_back(compute_ifft_row, std::ref(img), width, height, start, end, num_sub_threads);
    }
    compute_ifft_row(img, width, height, (num_threads - 1) * block_size, height, num_sub_threads);

    for (auto& th : threads) {
        th.join();
    }
}

void fft_image(
        std::vector<complex> &x, int width, int height,
        size_t num_sub_threads = 1
) {
    compute_fft_row(x, width, height, 0, height, num_sub_threads);
    compute_fft_column(x, width, height, 0, width, num_sub_threads);
}

void inverse_fft_image(
        std::vector<complex> &x, int width, int height,
        size_t num_sub_threads = 1
) {
    compute_ifft_row(x, width, height, 0, width, num_sub_threads);
    compute_ifft_column(x, width, height, 0, height, num_sub_threads);
}

// Function to compute the threshold value for compression
double get_threshold_value(const std::vector<complex>& data, double percentage_to_keep) {
    std::vector<double> abs_data;
    for (const auto& d : data) {
        double abs_value = std::abs(d);
        if (abs_value != 0.0) {  // Skip zero values
            abs_data.push_back(abs_value);
        }
    }
    if (abs_data.empty()) {
        std::cout << "No data to compress!" << std::endl;
        return 0.0;
    }
    std::sort(abs_data.begin(), abs_data.end());
    int index = static_cast<int>((1.0 - percentage_to_keep) * abs_data.size());
    index = std::max(0, std::min(index, static_cast<int>(abs_data.size() - 1)));
    return abs_data[index];
}

// Function to compress data using multiple threads
void compress_block_threads(std::vector<complex>& block, double cut_off, int start, int end) {
    for (int i = start; i < end; ++i) {
        if (std::abs(block[i]) < cut_off) {
            block[i] = 0;
        }
    }
}

void compress(std::vector<complex>& data, double cut_off, size_t num_threads) {
    int N = data.size();
    int block_size = N / num_threads;
    std::vector<std::thread> threads;

    for (int i = 0; i < num_threads - 1; ++i) {
        int start = i * block_size;
        int end = (i + 1) * block_size;
        threads.emplace_back(compress_block_threads, std::ref(data), cut_off, start, end);
    }
    compress_block_threads(data, cut_off, (num_threads - 1) * block_size, N);

    for (auto& th : threads) {
        th.join();
    }
}
//
//void compute_fft_row(std::vector<complex> &img, int width, int height, int start, int end, int num_sub_threads = 1) {
//    for (int i = start; i < end; i++) {
//        std::vector<complex> sliced(width);
//        for (int j = 0; j < width; j++) {
//            sliced[j] = img[i * width + j];
//        }
//
//        if (num_sub_threads == 1){
//            Radix2FFT(sliced);
//        }
//        else{
//            FFT_Radix2_Parallel(sliced, num_sub_threads);
//        }
//        for (int k = 0; k < sliced.size(); k++) {
//            img[i * width + k] = sliced[k];
//        }
//    }
//}
//
//void compute_fft_column(std::vector<complex> &img, int width, int height, int start, int end, int num_sub_threads = 1) {
//    for (int j = start; j < end; j++) {
//        std::vector<complex> sliced(height);
//        for (int i = 0; i < height; i++) {
//            sliced[i] = img[i * width + j];
//        }
//        if (num_sub_threads == 1){
//            Radix2FFT(sliced);
//        }
//        else{
//            FFT_Radix2_Parallel(sliced, num_sub_threads);
//        }
//        for (size_t k = 0; k < sliced.size(); k++) {
//            img[j + k * width] = sliced[k];
//        }
//    }
//}
//
//void compute_ifft_row(std::vector<complex> &img, int width, int height, int start, int end, int num_sub_threads = 1) {
//    for (int i = start; i < end; i++) {
//        std::vector<complex> sliced(width);
//        for (int j = 0; j < width; j++) {
//            sliced[j] = img[i * width + j];
//        }
//        if (num_sub_threads == 1){
//            InverseRadix2FFT(sliced);
//        }
//        else{
//            Inverse_Radix2_Parallel(sliced, num_sub_threads);
//        }
//        for (size_t k = 0; k < sliced.size(); k++) {
//            img[i * width + k] = sliced[k];
//        }
//    }
//}
//
//void compute_ifft_column(std::vector<complex> &img, int width, int height, int start, int end, int num_sub_threads = 1) {
//    for (int j = start; j < end; j++) {
//        std::vector<complex> sliced(height);
//        for (int i = 0; i < height; i++) {
//            sliced[i] = img[i * width + j];
//        }
//        if (num_sub_threads == 1){
//            InverseRadix2FFT(sliced);
//        }
//        else{
//            Inverse_Radix2_Parallel(sliced, num_sub_threads);
//        }
//        for (size_t k = 0; k < sliced.size(); k++) {
//            img[j + k * width] = sliced[k];
//        }
//    }
//}
//
//void parallel_fft_image(std::vector<complex> &img, int width, int height, int num_threads, bool double_parallel = false) {
//    int num_sub_threads = 1;
//    if (double_parallel) {
//        num_threads /= 2;
//        num_sub_threads = 2;
//    }
//    int block_size = height / num_threads;
//    std::vector<std::thread> threads(num_threads - 1);
//    int start = 0;
//    int end = block_size;
//    for (int i = 0; i < (num_threads - 1); i++) {
//        end = start+block_size;
//        threads[i] = std::thread(compute_fft_row,ref(img), width, height,start, end, num_sub_threads);
//        start = end;
//
//    }
//    compute_fft_row(img, width, height, start, height);
//    for (auto &th : threads) {
//        th.join();
//    }
//    block_size = width / num_threads;
//    start = 0;
//    for (int i = 0; i < (num_threads - 1); i++) {
//        end = start+block_size;
//        threads[i] = std::thread(compute_fft_column, ref(img), width, height, start, end, num_sub_threads);
//       start = end;
//    }
//    compute_fft_column(img, width, height, start, width);
//    for (auto &th : threads) {
//        th.join();
//    }
//}
//
//void parallel_inverse_fft_image(std::vector<complex> &x, int width, int height, int num_threads, bool double_parallel = false) {
//    int num_sub_threads = 1;
//    if (double_parallel) {
//        num_threads /= 2;
//        num_sub_threads = 2;
//    }
//    int block_size = width / num_threads;
//    std::vector<std::thread> threads(num_threads - 1);
//    int start = 0;
//    int end = block_size;
//    for (int i = 0; i < (num_threads - 1); i++) {
//        end = start+block_size;
//        threads[i] = std::thread(compute_ifft_column,ref(x), width, height,start, end, num_sub_threads);
//        start = end;
//    }
//    compute_ifft_column(x, width, height, start, width);
//    for (auto &th : threads) {
//        th.join();
//    }
//    block_size = height / num_threads;
//    start = 0;
//    for (int i = 0; i < (num_threads - 1); i++) {
//        end = start+block_size;
//        threads[i] = std::thread(compute_ifft_row,ref(x), width, height,start, end, num_sub_threads);
//        start = end;
//    }
//    compute_ifft_row(x, width, height, start, height);
//    for (auto &th : threads) {
//        th.join();
//    }
//}
//
//
//
//
//double get_threshold_value(std::vector<std::complex<double>>& data, double percentage_to_keep) {
//    std::vector<double> abs_data;
//    for (auto &d : data) {
//        double abs_value = abs(d);
//        if (abs_value != 0.0) {  // Skip zero values (in particular because of Radix2)
//            abs_data.push_back(abs_value);
//        }
//    }
//    if (abs_data.empty()) {
//        std::cout<< "No data to compress!" << std::endl;
//        return 0.0;
//    }
//    std::sort(abs_data.begin(), abs_data.end());
//    int index = static_cast<int>((1.0 - percentage_to_keep) * (abs_data.size() ));
//    index = std::max(0, std::min(index, static_cast<int>(abs_data.size())));
//    return abs_data[index];
//}
//
//void compress_block_threads(
//        std::vector<complex> &block, double cut_off, int start, int end
//) {
//    for (int i = start; i < end; i++) {
//        if (abs(block[i]) < cut_off){
//            block[i] = 0;
//        }
//    }
//}
//
//void compress(std::vector<complex> &data, double cut_off, size_t num_threads) {
//
//    int N = data.size();
//    int block_size = N/num_threads;
//    std::vector<std::thread> threads(num_threads - 1);
//    int start = 0;
//    for (int i = 0; i < (num_threads - 1); i++) {
//        int end = start + block_size;
//        threads[i] = std::thread(compress_block_threads,std::ref(data), cut_off, start, end);
//        start += block_size;
//    }
//    compress_block_threads(data, cut_off, start, N);
//    for(auto &th : threads) {
//        th.join();
//    }
//}