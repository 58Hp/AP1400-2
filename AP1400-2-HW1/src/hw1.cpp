#include <hw1.h>
#include <vector>
#include <stdexcept>
#include <random>
#include <iomanip>
#include <iostream>
#include <cmath>

using std::logic_error;
using std::vector;
namespace algebra {
	Matrix zeros(size_t n, size_t m) {
		Matrix matrix(n, vector<double>(m, 0));
		return matrix;
	}

	Matrix ones(size_t n, size_t m) {
		Matrix matrix(n, vector<double>(m, 1));
		return matrix;
	}

	Matrix random(size_t n, size_t m, double min, double max) {
		if (min > max) throw logic_error("�������Сֵ���ô������ֵ");
		Matrix matrix(n, vector<double>(m));
		std::mt19937 generator;
		std::uniform_real_distribution<double> distribution(min, max);
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				matrix[i][j] = distribution(generator);
			}
		}
		return matrix;
	}

	void show(const Matrix& matrix) {
		if (matrix.empty()) return;
		size_t n = matrix.size();
		size_t m = matrix[0].size();
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				std::cout << std::setprecision(3) << std::setw(6) << std::left << matrix[i][j];
			}
			std::cout << std::endl;
		}
	}

	Matrix multiply(const Matrix& matrix, double c) {
		if (matrix.empty()) return matrix;
		size_t n = matrix.size();
		size_t m = matrix[0].size();
		Matrix outMatrix(n, vector<double>(m));
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				outMatrix[i][j] = matrix[i][j] * c;
			}
		}
		return outMatrix;
	}

	Matrix multiply(const Matrix& matrix1, const Matrix& matrix2) {
		if (matrix1.empty()) return {};
		if (matrix2.empty()) return {};
		size_t row1 = matrix1.size();
		size_t col1 = matrix1[0].size();
		size_t row2 = matrix2.size();
		size_t col2 = matrix2[0].size();
		if (col1 != row2) throw logic_error("������״��ƥ�䣡");
		Matrix outMatrix(row1, vector<double>(col2));
		for (int i = 0; i < row1; i++) {
			for (int j = 0; j < col2; j++) {
				double temp = 0;
				for (int k = 0; k < col1; k++) {
					temp += matrix1[i][k] * matrix2[k][j];
				}
				outMatrix[i][j] = temp;
			}
		}
		return outMatrix;
	}

	Matrix sum(const Matrix& matrix, double c) {
		if (matrix.empty()) return {};
		size_t n = matrix.size();
		size_t m = matrix[0].size();
		Matrix outMatrix(n, vector<double>(m));
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				outMatrix[i][j] = matrix[i][j] + c;
			}
		}
		return outMatrix;
	}

	Matrix sum(const Matrix& matrix1, const Matrix& matrix2) {
		if (matrix1.empty() && matrix2.empty()) return {};
		if (matrix1.empty() || matrix2.empty()) throw logic_error("����ͬ�ͣ�");
		size_t row1 = matrix1.size();
		size_t col1 = matrix1[0].size();
		size_t row2 = matrix2.size();
		size_t col2 = matrix2[0].size();
		if (row1 != row2 || col1 != col2) throw logic_error("����ͬ�ͣ�");
		Matrix outMatrix(row1, vector<double>(col1));
		for (int i = 0; i < row1; i++) {
			for (int j = 0; j < col2; j++) {
				outMatrix[i][j] = matrix1[i][j] + matrix2[i][j];
			}
		}
		return outMatrix;
	}

	Matrix transpose(const Matrix& matrix) {
		if (matrix.empty()) return {};
		size_t row = matrix.size();
		size_t col = matrix[0].size();
		Matrix outMatrix(col, vector<double>(row));
		for (int i = 0; i < row; i++) {
			for (int j = 0; j < col; j++) {
				outMatrix[j][i] = matrix[i][j];
			}
		}
		return outMatrix;
	}

	Matrix minor(const Matrix& matrix, size_t n, size_t m) {
		if (matrix.empty()) return {};
		size_t row = matrix.size();
		size_t col = matrix[0].size();
		if (n < 0 || m < 0 || n >= row || m >= col) throw logic_error("����ʽ�Ƿ���");
		Matrix outMatrix(row - 1, vector<double>(col - 1));
		if (row - 1 <= 0 || col - 1 <= 0) {
			return outMatrix;
		}
		for (int i = 0; i < row; i++) {
			for (int j = 0; j < col; j++) {
				if (i == n || j == m) {
					continue;
				}
				else if (i > n && j > m) {
					outMatrix[i - 1][j - 1] = matrix[i][j];
				}
				else if (i < n && j > m) {
					outMatrix[i][j - 1] = matrix[i][j];
				}
				else if (i > n && j < m) {
					outMatrix[i - 1][j] = matrix[i][j];
				}
				else if (i < n && j < m) {
					outMatrix[i][j] = matrix[i][j];
				}
			}
		}
		return outMatrix;
	}

	double determinant(const Matrix& matrix) {
		if (matrix.empty()) return 1;//�վ��������ʽԼ��Ϊ1
		size_t n = matrix.size();
		size_t m = matrix[0].size();
		if (m != n) throw logic_error("�����Ƿ�������޷���������ʽ��");
		if (n == 1) {
			return matrix[0][0];
		}
		else if (n == 2) {
			return matrix[0][0] * matrix[1][1] - matrix[1][0] * matrix[0][1];
		}
		double sum = 0;
		for (int i = 0; i < n; i++) {
			sum += matrix[i][0] * pow(-1, i) * determinant(minor(matrix, i, 0));
		}
		return sum;
	}

	Matrix inverse(const Matrix& matrix) {
		if (matrix.empty()) return {};//�վ���������Ϊ�վ���
		size_t n = matrix.size();
		size_t m = matrix[0].size();
		if (m != n) {
			throw logic_error("�����Ƿ�������޷����������");
		}
		else if (determinant(matrix) == 0) {
			throw logic_error("�þ��󲻴��������");
		}
		double deter = determinant(matrix);
		Matrix outMatrix(m, std::vector<double>(n));
		for (int i = 0; i < n; i++) {//���������
			for (int j = 0; j < m; j++) {
				outMatrix[j][i] = pow(-1, i + j) * determinant(minor(matrix, i, j));
				outMatrix[j][i] /= deter;
			}
		}
		return outMatrix;
	}

	Matrix concatenate(const Matrix& matrix1, const Matrix& matrix2, int axis) {
		if (matrix1.empty() && matrix2.empty()) return {};
		if (matrix1.empty() || matrix2.empty()) throw logic_error("���ڿվ����µ��޷����ӣ�");
		size_t row1 = matrix1.size();
		size_t col1 = matrix1[0].size();
		size_t row2 = matrix2.size();
		size_t col2 = matrix2[0].size();
		if (axis == 0) {
			if (col1 != col2) throw logic_error("�����Ȳ�ƥ�䣬�޷������������ӣ�");
			Matrix outMatrix(row1 + row2, vector<double>(col1));
			for (int i = 0; i < row1; i++) {
				for (int j = 0; j < col1; j++) {
					outMatrix[i][j] = matrix1[i][j];
				}
			}
			for (int i = 0; i < row2; i++) {
				for (int j = 0; j < col2; j++) {
					outMatrix[row1 + i][j] = matrix2[i][j];
				}
			}
			return outMatrix;
		}
		else if (axis == 1) {
			if (row1 != row2) throw logic_error("����߶Ȳ�ƥ�䣬�޷����к������ӣ�");
			Matrix outMatrix(row1, vector<double>(col1 + col2));
			for (int i = 0; i < row1; i++) {
				for (int j = 0; j < col1; j++) {
					outMatrix[i][j] = matrix1[i][j];
				}
			}
			for (int i = 0; i < row2; i++) {
				for (int j = 0; j < col2; j++) {
					outMatrix[i][col1 + j] = matrix2[i][j];
				}
			}
			return outMatrix;
		}
		else throw logic_error("δ֪axis������");
	}

	Matrix ero_swap(const Matrix& matrix, size_t r1, size_t r2) {
		if (matrix.empty()) throw logic_error("�޷��Կվ�������в�����");
		size_t row = matrix.size();
		size_t col = matrix[0].size();
		if (r1 < 0 || r1 >= row || r2 < 0 || r2 >= row) throw logic_error("���󲻰������У�");
		Matrix outMatrix{ matrix };
		outMatrix[r1].swap(outMatrix[r2]);
		return outMatrix;
	}

	Matrix ero_multiply(const Matrix& matrix, size_t r, double c) {
		if (matrix.empty()) return {};
		size_t row = matrix.size();
		size_t col = matrix[0].size();
		if (r < 0 || r >= row) throw logic_error("���󲻰������У�");
		Matrix outMatrix{ matrix };
		for (int i = 0; i < row; i++) {
			outMatrix[r][i] *= c;
		}
		return outMatrix;
	}

	Matrix ero_sum(const Matrix& matrix, size_t r1, double c, size_t r2) {
		if (matrix.empty()) throw logic_error("�޷��Կվ�������в�����");
		size_t row = matrix.size();
		size_t col = matrix[0].size();
		if (r1 < 0 || r1 >= row || r2 < 0 || r2 >= row) throw logic_error("���󲻰������У�");
		Matrix outMatrix{ matrix };
		for (int i = 0; i < col; i++) {
			outMatrix[r2][i] += outMatrix[r1][i] * c;
		}
		return outMatrix;
	}

	Matrix upper_triangular(const Matrix& matrix) {
		if (matrix.empty()) return {};
		size_t row = matrix.size();
		size_t col = matrix[0].size();
		if (row != col) throw logic_error("�˾��󲻴��������Ǿ���");//����Ҫ��
		Matrix outMatrix{ matrix };
		for (int i = 0; i < row && i < col; i++) {//ÿ����һ��Ϊ��λ����
			if (outMatrix[i][i] == 0) {//����Ԫ��Ϊ0��������Ԫ�ط����н��������ȫΪ0�򲻴���
				int j;
				//������Ԫ�ط����н���
				for (j = i + 1; j < row; j++) {
					if (outMatrix[j][i] != 0) break;
				}
				if (j == row) continue;//δ�ҵ��ɽ�������
				outMatrix = ero_swap(outMatrix, i, j);
			}
			for (int j = i + 1; j < row; j++) {
				outMatrix = ero_sum(outMatrix, i, -outMatrix[j][i] / outMatrix[i][i], j);
			}
		}
		return outMatrix;
	}
}