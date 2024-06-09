#include "plot.cpp"
#define STB_IMAGE_IMPLEMENTATION
#include "STB/stb_image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "STB/stb_image_write.h"


/*------------------------ Image read ------------------------*/

std::vector<complex>* read_img(std::string& filename, int &width, int &height) {
    int channels;
    uint8_t* rgb_img = stbi_load(filename.c_str(), &width, &height, &channels, 3);
    if (rgb_img == nullptr) {
        std::cerr << "Error loading image: " << filename << std::endl;
        return nullptr;
    }

    std::cout << "Image Size: " << width << "x" << height <<
              " (" << channels << " channels)" << std::endl;

    size_t rgb_img_size = width * height * 3; // We are loading the image with 3 channels (RGB)
    size_t grey_img_size = width * height;

    std::vector<complex> img_data(grey_img_size);

    size_t k = 0;
    for (uint8_t *p = rgb_img; p < rgb_img + rgb_img_size; p += 3) { // Move by 3 for RGB
        img_data[k] = 0.299 * (*p) + 0.587 * (*(p + 1)) + 0.114 * (*(p + 2)); // Grayscale conversion
        k++;
    }

    stbi_image_free(rgb_img);

    std::cout << "Image read successfully!" << std::endl;
    return new std::vector<complex>(std::move(img_data));
}

/*------------------------ Image write ------------------------*/
void write_img(std::vector<complex> &img_data, int width, int height, std::string filename) {
    std::cout << "Writing image to: " << filename << std::endl;

    size_t out_img_size = width * height;
    std::vector<uint8_t> out_img(out_img_size);

    for (size_t k = 0; k < out_img_size; k++) {
        out_img[k] = static_cast<uint8_t>(std::max(std::min(255.0, img_data[k].real()), 0.0));
    }

    stbi_write_png(filename.c_str(), width, height, 1, out_img.data(), width * sizeof(uint8_t));
}


/*------------------------ Image FFT ------------------------*/



void compute_fft_row(std::vector<complex> &img, int width, int height, int start, int end, int num_sub_threads = 1) {
    for (int i = start; i < end; i++) {
        std::vector<complex> sliced(width);
        for (int j = 0; j < width; j++) {
            sliced[j] = img[i * width + j];
        }

        if (num_sub_threads == 1){
            Radix2FFT(sliced);
        }
        else{
            FFT_Radix2_Parallel(sliced, num_sub_threads);
        }
        for (int k = 0; k < sliced.size(); k++) {
            img[i * width + k] = sliced[k];
        }
    }
}

void compute_fft_column(std::vector<complex> &img, int width, int height, int start, int end, int num_sub_threads = 1) {
    for (int j = start; j < end; j++) {
        std::vector<complex> sliced(height);
        for (int i = 0; i < height; i++) {
            sliced[i] = img[i * width + j];
        }
        if (num_sub_threads == 1){
            Radix2FFT(sliced);
        }
        else{
            FFT_Radix2_Parallel(sliced, num_sub_threads);
        }
        for (size_t k = 0; k < sliced.size(); k++) {
            img[j + k * width] = sliced[k];
        }
    }
}

void compute_ifft_row(std::vector<complex> &img, int width, int height, int start, int end, int num_sub_threads = 1) {
    for (int i = start; i < end; i++) {
        std::vector<complex> sliced(width);
        for (int j = 0; j < width; j++) {
            sliced[j] = img[i * width + j];
        }
        if (num_sub_threads == 1){
            InverseRadix2FFT(sliced);
        }
        else{
            Inverse_Radix2_Parallel(sliced, num_sub_threads);
        }
        for (size_t k = 0; k < sliced.size(); k++) {
            img[i * width + k] = sliced[k];
        }
    }
}

void compute_ifft_column(std::vector<complex> &img, int width, int height, int start, int end, int num_sub_threads = 1) {
    for (int j = start; j < end; j++) {
        std::vector<complex> sliced(height);
        for (int i = 0; i < height; i++) {
            sliced[i] = img[i * width + j];
        }
        if (num_sub_threads == 1){
            InverseRadix2FFT(sliced);
        }
        else{
            Inverse_Radix2_Parallel(sliced, num_sub_threads);
        }
        for (size_t k = 0; k < sliced.size(); k++) {
            img[j + k * width] = sliced[k];
        }
    }
}

void parallel_fft_image(std::vector<complex> &img, int width, int height, int num_threads, bool double_parallel = false) {
    int num_sub_threads = 1;
    if (double_parallel) {
        num_threads /= 2;
        num_sub_threads = 2;
    }
    int block_size = height / num_threads;
    std::vector<std::thread> threads(num_threads - 1);
    int start = 0;
    int end = block_size;
    for (int i = 0; i < (num_threads - 1); i++) {
        end = start+block_size;
        threads[i] = std::thread(compute_fft_row,ref(img), width, height,start, end, num_sub_threads);
        start = end;

    }
    compute_fft_row(img, width, height, start, height);
    for (auto &th : threads) {
        th.join();
    }
    block_size = width / num_threads;
    start = 0;
    for (int i = 0; i < (num_threads - 1); i++) {
        end = start+block_size;
        threads[i] = std::thread(compute_fft_column, ref(img), width, height, start, end, num_sub_threads);
       start = end;
    }
    compute_fft_column(img, width, height, start, width);
    for (auto &th : threads) {
        th.join();
    }
}

void parallel_inverse_fft_image(std::vector<complex> &x, int width, int height, int num_threads, bool double_parallel = false) {
    int num_sub_threads = 1;
    if (double_parallel) {
        num_threads /= 2;
        num_sub_threads = 2;
    }
    int block_size = width / num_threads;
    std::vector<std::thread> threads(num_threads - 1);
    int start = 0;
    int end = block_size;
    for (int i = 0; i < (num_threads - 1); i++) {
        end = start+block_size;
        threads[i] = std::thread(compute_ifft_column,ref(x), width, height,start, end, num_sub_threads);
        start = end;
    }
    compute_ifft_column(x, width, height, start, width);
    for (auto &th : threads) {
        th.join();
    }
    block_size = height / num_threads;
    start = 0;
    for (int i = 0; i < (num_threads - 1); i++) {
        end = start+block_size;
        threads[i] = std::thread(compute_ifft_row,ref(x), width, height,start, end, num_sub_threads);
        start = end;
    }
    compute_ifft_row(x, width, height, start, height);
    for (auto &th : threads) {
        th.join();
    }
}




double get_threshold_value(std::vector<std::complex<double>>& data, double percentage_to_keep) {
    std::vector<double> abs_data;
    for (auto &d : data) {
        double abs_value = abs(d);
        if (abs_value != 0.0) {  // Skip zero values (in particular because of Radix2)
            abs_data.push_back(abs_value);
        }
    }
    if (abs_data.empty()) {
        std::cout<< "No data to compress!" << std::endl;
        return 0.0;
    }
    std::sort(abs_data.begin(), abs_data.end());
    int index = static_cast<int>((1.0 - percentage_to_keep) * (abs_data.size() ));
    index = std::max(0, std::min(index, static_cast<int>(abs_data.size())));
    return abs_data[index];
}

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