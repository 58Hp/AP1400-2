#include <hw1.h>
#include <vector>
#include <stdexcept>
#include <random>
#include <iomanip>
#include <iostream>
#include <cmath>

Matrix algebra::zeros(size_t n, size_t m) {
	if (n == 0 || m == 0) {
		throw std::logic_error("矩阵长宽不得等于0！");
	}
	Matrix matrix(n, std::vector<double>(m,0));
	return matrix;
}

Matrix algebra::ones(size_t n, size_t m) {
	if (n == 0 || m == 0) {
		throw std::logic_error("矩阵长宽不得等于0！");
	}
	Matrix matrix(n, std::vector<double>(m, 1));
	return matrix;
}

Matrix algebra::random(size_t n, size_t m, double min, double max) {
	if (n == 0 || m == 0) {
		throw std::logic_error("矩阵长宽不得等于0！");
	}
	if (min > max) {
		throw std::logic_error("随机数最小值不得大于最大值");
	}
	Matrix matrix(n, std::vector<double>(m));
	std::mt19937 generator;
	std::uniform_real_distribution<double> distribution(min, max);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			matrix[i][j] = distribution(generator);
		}
	}
	return matrix;
}

void algebra::show(const Matrix& matrix) {
	size_t n = matrix.size();
	if (n <= 0) {
		return;
	}
	size_t m = matrix[0].size();
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			std::cout << std::setprecision(3) << std::setw(6)<< std::left << matrix[i][j];
		}
		std::cout << std::endl;
	}
}

Matrix algebra::multiply(const Matrix& matrix, double c) {
	size_t n = matrix.size();
	if (n <= 0) {
		return matrix;
	}
	size_t m = matrix[0].size();
	Matrix outMatrix(n, std::vector<double>(m));
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			outMatrix[i][j] = matrix[i][j] * c;
		}
	}
	return outMatrix;
}

Matrix algebra::multiply(const Matrix& matrix1, const Matrix& matrix2) {
	size_t n1 = matrix1.size();
	size_t n2 = matrix2.size();
	if (n1 <= 0) {
		return matrix1;
	}
	else if(n2 <= 0){
		return matrix2;
	}
	size_t m1 = matrix1[0].size();
	size_t m2 = matrix2[0].size();
	if (m1 != n2) {
		throw std::logic_error("矩阵形状不匹配！");
	}
	Matrix outMatrix(n1, std::vector<double>(m2));
	for (int i = 0; i < n1; i++) {
		for (int j = 0; j < m2; j++) {
			double temp = 0;
			for (int n = 0; n < m1; n++) {
				temp += matrix1[i][n] * matrix2[n][j];
			}
			outMatrix[i][j] = temp;
		}
	}
	return outMatrix;
}

Matrix  algebra::sum(const Matrix& matrix, double c) {
	size_t n = matrix.size();
	if (n <= 0) {
		return matrix;
	}
	size_t m = matrix[0].size();
	Matrix outMatrix(n, std::vector<double>(m));
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			outMatrix[i][j] = matrix[i][j] + c;
		}
	}
	return outMatrix;
}

Matrix algebra::sum(const Matrix& matrix1, const Matrix& matrix2) {
	size_t n1 = matrix1.size();
	size_t n2 = matrix2.size();
	if (n1 != n2) {
		throw std::logic_error("矩阵不同型！");
	}
	if (n1 <= 0) {
		return matrix1;
	}
	else if (n2 <= 0) {
		return matrix2;
	}
	size_t m1 = matrix1[0].size();
	size_t m2 = matrix2[0].size();
	if (m1 != m2) {
		throw std::logic_error("矩阵不同型！");
	}
	Matrix outMatrix(n1, std::vector<double>(m1));
	for (int i = 0; i < n1; i++) {
		for (int j = 0; j < n2; j++) {
			outMatrix[i][j] = matrix1[i][j]+ matrix2[i][j];
		}
	}
	return outMatrix;
}

Matrix algebra::transpose(const Matrix& matrix) {
	size_t n = matrix.size();
	if (n <= 0) {
		return matrix;
	}
	size_t m = matrix[0].size();
	Matrix outMatrix(m, std::vector<double>(n));
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			outMatrix[j][i] = matrix[i][j];
		}
	}
	return outMatrix;
}

Matrix algebra::minor(const Matrix& matrix, size_t n, size_t m) {
	size_t height = matrix.size();
	if (height <= 0) {
		throw std::logic_error("矩阵不能为空！");
	}
	size_t width = matrix[0].size();
	if (n < 0 || m< 0 || n >= height || m>= width) {
		throw std::logic_error("余子式非法！");
	}
	Matrix outMatrix(height-1, std::vector<double>(width-1));
	if (height - 1 <= 0 || width - 1 <= 0) {
		return outMatrix;
	}
	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			if (i == n || j == m) {
				continue;
			}
			else if (i > n && j > m) {
				outMatrix[i-1][j-1] = matrix[i][j];
			}
			else if (i < n && j > m) {
				outMatrix[i][j-1] = matrix[i][j];
			}
			else if (i > n && j < m) {
				outMatrix[i-1][j] = matrix[i][j];
			}
			else if (i < n && j < m) {
				outMatrix[i][j] = matrix[i][j];
			}
		}
	}
	return outMatrix;
}

double algebra::determinant(const Matrix& matrix) {
	size_t n = matrix.size();
	if (n <= 0) {
		return 1;
	}
	size_t m = matrix[0].size();
	if (m != n) {
		throw std::logic_error("无法计算行列式！");
	}
	if (n == 1) {
		return matrix[0][0];
	}
	else if (n == 2) {
		return matrix[0][0]* matrix[1][1] - matrix[1][0] * matrix[0][1];
	}
	double sum=0;
	for (int i = 0; i < n; i++) {
		sum +=matrix[i][0]*pow(-1,i)*algebra::determinant(algebra::minor(matrix, i, 0));
	}
	return sum;
}

Matrix algebra::inverse(const Matrix& matrix) {
	size_t n = matrix.size();
	if (n <= 0) {//空矩阵的逆矩阵为空矩阵
		return matrix;
	}
	size_t m = matrix[0].size();
	if (m != n) {
		throw std::logic_error("无法计算逆矩阵！");
	}
	if (algebra::determinant(matrix) == 0) {
		throw std::logic_error("该矩阵不存在逆矩阵");
	}
	double deter = algebra::determinant(matrix);
	Matrix outMatrix(m, std::vector<double>(n));
	for (int i = 0; i < n; i++) {//计算逆矩阵
		for (int j = 0; j < m; j++) {
			outMatrix[j][i] = pow(-1,i+j)*algebra::determinant(algebra::minor(matrix, i, j));
			outMatrix[j][i] /= deter;
		}
	}
	return outMatrix;
}