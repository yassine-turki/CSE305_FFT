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
void compute_fft_row(std::vector<complex>& img, int width, int height, int start, int end, size_t num_threads = 10) {
    for (int i = start; i < end; ++i) {
        std::vector<complex> row(img.begin() + i * width, img.begin() + (i + 1) * width);

        FFT_Radix2_Parallel(row, num_threads);

        std::copy(row.begin(), row.end(), img.begin() + i * width);
    }
}


void compute_fft_column(std::vector<complex>& img, int width, int height, int start, int end, size_t num_threads = 10) {
    for (int j = start; j < end; ++j) {
        std::vector<complex> column(height);
        for (int i = 0; i < height; ++i) {
            column[i] = img[i * width + j];
        }


        FFT_Radix2_Parallel(column, num_threads);


        for (int i = 0; i < height; ++i) {
            img[i * width + j] = column[i];
        }
    }
}

// Function to compute inverse FFT for rows
void compute_ifft_row(std::vector<complex>& img, int width, int height, int start, int end, size_t num_threads = 10) {
    for (int i = start; i < end; ++i) {
        std::vector<complex> row(img.begin() + i * width, img.begin() + (i + 1) * width);


        Inverse_Radix2_Parallel(row, num_threads);


        std::copy(row.begin(), row.end(), img.begin() + i * width);
    }
}

// Function to compute inverse FFT for columns
void compute_ifft_column(std::vector<complex>& img, int width, int height, int start, int end, size_t num_threads = 10) {
    for (int j = start; j < end; ++j) {
        std::vector<complex> column(height);
        for (int i = 0; i < height; ++i) {
            column[i] = img[i * width + j];
        }
        Inverse_Radix2_Parallel(column, num_threads);
        for (int i = 0; i < height; ++i) {
            img[i * width + j] = column[i];
        }
    }
}

void compute_fft_img(std::vector<complex> &x, int width, int height, size_t num_threads = 10) {
    compute_fft_row(x, width, height, 0, height, num_threads);
    compute_fft_column(x, width, height, 0, width, num_threads);
}

void compute_ifft_img(std::vector<complex> &x, int width, int height, size_t num_threads = 10) {
    compute_ifft_row(x, width, height, 0, width, num_threads);
    compute_ifft_column(x, width, height, 0, height, num_threads);
}


/*------------------------ Image Compression ------------------------*/

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
void compress_block_threads(std::vector<complex>& block, double threshold, int start, int end) {
    for (int i = start; i < end; ++i) {
        if (std::abs(block[i]) < threshold) {
            block[i] = 0.;
        }
    }
}

void compress(std::vector<complex>& data, double threshold, size_t num_threads) {
    int N = data.size();
    int block_size = N / num_threads;
    std::vector<std::thread> threads;

    for (int i = 0; i < num_threads - 1; ++i) {
        int start = i * block_size;
        int end = (i + 1) * block_size;
        threads.emplace_back(compress_block_threads, std::ref(data), threshold, start, end);
    }
    compress_block_threads(data, threshold, (num_threads - 1) * block_size, N);

    for (auto& th : threads) {
        th.join();
    }
}
