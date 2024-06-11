# CSE305_FFT

## Overview
This project explores the implementation and optimization of the Fast Fourier Transform (FFT) using the Radix-2 algorithm. The project focuses on parallelizing the different steps involved in computing the FFT, as well as auxiliary functions, to enhance efficiency. Additionally, the project applies the FFT algorithm to various real-world use cases such as image compression, signal compression, and polynomial multiplication.

## Repository Contents
The repository consists of several C++ files, a Jupyter notebook for plotting benchmark data, and datasets and images used in the project.

### C++ Files
- **fft.cpp:** Contains implementations for Naive DFT and IDFT, their parallel versions, and the Radix-2 algorithm (FFT and IFFT).
- **poly_ops.cpp:** Includes operations for polynomial manipulation, including naive and FFT-based polynomial multiplication, along with their parallel implementations.
- **ntt.cpp:** Contains operations for polynomials in the quotient ring, including the Number Theoretic Transform (NTT).
- **data.cpp:** Functions for generating synthetic data, reading data from CSV files, and writing data to CSV files.
- **benchmark.cpp:** Functions for benchmarking the computation time of the FFT implementations.
- **compress.cpp:** Functions for reading and writing images, computing the FFT for images, and compressing data.
- **test.cpp**: Various testing functions for validating the implemented algorithms and compression techniques.
- **main.cpp** Main file with booleans to run different parts of the code (polynomial multiplication, FFT on weather data, image compression, etc.).
### Additional Contents
- STBI Folder: Handles reading and writing images (adapted from the computer graphics class CSE306).
- Plot_benchmark_fft-2.ipynb: Jupyter notebook for plotting benchmark data.
- data Folder: Contains the weather dataset used for the project.
- images Folder: Contains the images used for testing image compression.
## Applications
### Benchmarking
We benchmarked the performance of the FFT algorithms using a small and large weather datasets. The results demonstrate the efficiency of the Radix-2 algorithm compared to the naive DFT, especially with parallel implementations.

### Compression
#### Signal Compression
We applied the FFT to compress weather data by retaining only significant frequency components. The results show that significant data reduction is possible while preserving the essential features of the signal.

#### Image Compression
We used the FFT to compress images by converting them to the frequency domain, filetring negligible frequencies, and reconstructing the image using the inverse FFT. The results illustrate a great compression quality with minimal loss of detail.

#### Polynomial Multiplication
We implemented FFT-based and NTT-based polynomial multiplication algorithms. The FFT approach is used for polynomials with complex coefficients, while the NTT is used for polynomials with integer coefficients modulo a prime number.

## Conclusion
This project demonstrates the efficiency and applicability of the FFT and its parallel implementation in various real-world scenarios. 

## Authors:
Alexandru Serban & Yassine Turki
