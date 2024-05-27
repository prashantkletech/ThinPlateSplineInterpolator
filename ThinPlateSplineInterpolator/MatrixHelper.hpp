#include <vector>
#include <cmath>
#include <stdexcept>

template <typename T>
std::vector<std::vector<T>> transposeMatrix(const std::vector<std::vector<T>>& matrix) {
    std::vector<std::vector<T>> transposed(matrix[0].size(), std::vector<T>(matrix.size()));
    for (size_t i = 0; i < matrix.size(); ++i) {
        for (size_t j = 0; j < matrix[0].size(); ++j) {
            transposed[j][i] = matrix[i][j];
        }
    }
    return transposed;
}

template <typename T>
std::vector<std::vector<T>> multiplyMatrices(const std::vector<std::vector<T>>& a, const std::vector<std::vector<T>>& b) {
    std::vector<std::vector<T>> result(a.size(), std::vector<T>(b[0].size()));
    for (size_t i = 0; i < a.size(); ++i) {
        for (size_t j = 0; j < b[0].size(); ++j) {
            result[i][j] = 0;
            for (size_t k = 0; k < b.size(); ++k) {
                result[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    return result;
}

template <typename T>
std::vector<T> multiplyMatrixVector(const std::vector<std::vector<T>>& matrix, const std::vector<T>& vec) {
    std::vector<T> result(matrix.size(), 0);
    for (size_t i = 0; i < matrix.size(); ++i) {
        for (size_t j = 0; j < vec.size(); ++j) {
            result[i] += matrix[i][j] * vec[j];
        }
    }
    return result;
}

template <typename T>
std::vector<std::vector<T>> invertMatrix(const std::vector<std::vector<T>>& matrix) {
    int n = matrix.size();
    std::vector<std::vector<T>> augmentedMatrix(n, std::vector<T>(2 * n));

    // Create the augmented matrix [A|I]
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            augmentedMatrix[i][j] = matrix[i][j];
        }
        augmentedMatrix[i][n + i] = 1;
    }

    // Perform Gaussian elimination
    for (int i = 0; i < n; ++i) {
        // Find the pivot row and swap
        int maxRow = i;
        for (int k = i + 1; k < n; ++k) {
            if (std::abs(augmentedMatrix[k][i]) > std::abs(augmentedMatrix[maxRow][i])) {
                maxRow = k;
            }
        }
        std::swap(augmentedMatrix[i], augmentedMatrix[maxRow]);

        // Make sure the matrix is invertible
        if (std::abs(augmentedMatrix[i][i]) < 1e-10) {
            throw std::runtime_error("Matrix is singular and cannot be inverted.");
        }

        // Normalize the pivot row
        T pivot = augmentedMatrix[i][i];
        for (int j = 0; j < 2 * n; ++j) {
            augmentedMatrix[i][j] /= pivot;
        }

        // Eliminate the current column
        for (int k = 0; k < n; ++k) {
            if (k != i) {
                T factor = augmentedMatrix[k][i];
                for (int j = 0; j < 2 * n; ++j) {
                    augmentedMatrix[k][j] -= factor * augmentedMatrix[i][j];
                }
            }
        }
    }

    // Extract the inverse matrix from the augmented matrix
    std::vector<std::vector<T>> inverseMatrix(n, std::vector<T>(n));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            inverseMatrix[i][j] = augmentedMatrix[i][n + j];
        }
    }

    return inverseMatrix;
}
