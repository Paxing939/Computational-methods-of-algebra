#include <cmath>
#include <iomanip>
#include <utility>
#include "Functions.h"

class Matrix {

public:
	Matrix() {
		srand(time(NULL));
		// random size of matrix
		int n = rand() % 3 + 12;

		std::vector<int> no_diag;
		for (int i = 0; i < (n * n - n) / 2; i++) {
			no_diag.push_back(rand() % 5 - 4);
		}

		matrix_k_2_.resize(n);
		for (auto& line : matrix_k_2_) {
			line.resize(n, -100);
		}

		// fill matrix
		int p = 0;
		double sum_1_col = 0.0f, sum_1_line = 0.0f;
		for (int i = 0; i < matrix_k_2_.size(); i++) {
			for (int j = 0; j < matrix_k_2_.size(); j++) {
				if (i != j) {
					matrix_k_2_[i][j] = rand() % 5 - 4;
				}
				if (i == 0 && j != 0) {
					sum_1_line += matrix_k_2_[0][j];
				}
				if (j == 0 && i != 0) {
					sum_1_col += matrix_k_2_[i][0];
				}
			}
		}

		matrix_k_0_ = matrix_k_2_;

		matrix_k_0_[0][0] = -sum_1_line + 1;
		matrix_k_2_[0][0] = -sum_1_line + pow(10, -2);

		double sum_1_col_k_0 = 0.0f, sum_1_col_k_2 = 0.0f;
		for (int i = 1; i < matrix_k_0_.size(); i++) {
			sum_1_col_k_0 = 0.0f;
			sum_1_col_k_2 = 0.0f;
			for (int j = 0; j < matrix_k_0_.size(); j++) {
				if (i != j) {
					sum_1_col_k_0 += matrix_k_0_[i][j];
					sum_1_col_k_2 += matrix_k_2_[i][j];
				}
			}
			matrix_k_0_[i][i] = sum_1_col - matrix_k_0_[i][0];
		}

		int m = 17;
		for (int i = 0; i < n; i++) {
			x_.push_back(m + i);
		}

		// multiply x_ on matrix to get vector b
		// if after Gauss method answer vector and x_ will be the same
		// then we will know method works correctly
		b_k_0_ = MultMat(matrix_k_0_, x_);
		b_k_2_ = MultMat(matrix_k_2_, x_);

		for (int i = 0; i < n; i++) {
			matrix_k_0_[i].resize(n + 1);
			matrix_k_2_[i].resize(n + 1);
			matrix_k_0_[i][n] = b_k_0_[i];
			matrix_k_2_[i][n] = b_k_2_[i];
		}

		// write matrices in the logs.txt
		LogMatrix(matrix_k_0_);
		LogMatrix(matrix_k_2_);
	}


	std::pair<std::vector<double>, std::vector<double>> GaussWithoutMainElement() {
		auto matrix_k_0 = matrix_k_0_, matrix_k_2 = matrix_k_2_;
		int line_length = matrix_k_0.size();

		// forward Gauss method
		for (int i = 0; i < line_length; ++i) {
			for (int s = i + 1; s < line_length; ++s) {
				double h_k_0 = matrix_k_0[s][i] / matrix_k_0[i][i];
				double h_k_2 = matrix_k_2[s][i] / matrix_k_2[i][i];
				for (int j = i; j < line_length + 1; ++j) {
					matrix_k_0[s][j] -= matrix_k_0[i][j] * h_k_0;
					matrix_k_2[s][j] -= matrix_k_2[i][j] * h_k_2;
				}
			}
		}

		// back Gauss method
		std::vector<double> vals_k_0(line_length), vals_k_2(line_length);
		double x_k_0 = 0, x_k_2 = 0;
		for (int i = line_length - 1; i >= 0; --i) {
			if (i > line_length) {
				break;
			}

			x_k_0 = matrix_k_0[i][line_length];
			x_k_2 = matrix_k_2[i][line_length];
			for (int j = line_length - 1; j > i; --j) {
				if (i == line_length - 1 && j == line_length - 1) {
					vals_k_0[i] = matrix_k_0[i][j];
					vals_k_2[i] = matrix_k_2[i][j];
					continue;
				}
				else {
					x_k_0 -= matrix_k_0[i][j] * vals_k_0[j];
					x_k_2 -= matrix_k_2[i][j] * vals_k_2[j];
				}
			}
			vals_k_0[i] = x_k_0 / matrix_k_0[i][i];
			vals_k_2[i] = x_k_2 / matrix_k_2[i][i];
		}

		return { vals_k_0, vals_k_2 };
	}

	std::pair<std::vector<double>, std::vector<double>> GaussWithMainElement() {
		auto matrix_k_0 = matrix_k_0_, matrix_k_2 = matrix_k_2_;
		int line_length = matrix_k_0.size();

		// forward Gauss method with choosing of main element
		int max_index_k_0 = 0, max_index_k_2 = 0;
		double max_k_0 = 0, max_k_2 = 0;
		for (int k = 0; k < line_length; ++k) {
			max_index_k_0 = k;
			max_k_0 = 0;
			max_index_k_2 = k;
			max_k_2 = 0;

			for (int i = k; i < line_length; ++i) {
				if (matrix_k_0[i][k] > max_k_0) {
					max_k_0 = matrix_k_0[i][k];
					max_index_k_0 = i;
				}

				if (matrix_k_2[i][k] > max_k_2) {
					max_k_2 = matrix_k_2[i][k];
					max_index_k_2 = i;
				}
			}
			swap(matrix_k_0[k], matrix_k_0[max_index_k_0]);
			swap(matrix_k_2[k], matrix_k_2[max_index_k_2]);
		}

		for (int i = 0; i < line_length; ++i) {
			for (int s = i + 1; s < line_length; ++s) {
				double h_k_0 = matrix_k_0[s][i] / matrix_k_0[i][i];
				double h_k_2 = matrix_k_2[s][i] / matrix_k_2[i][i];

				for (int j = i; j < line_length + 1; ++j) {
					matrix_k_0[s][j] -= matrix_k_0[i][j] * h_k_0;
					matrix_k_2[s][j] -= matrix_k_2[i][j] * h_k_2;
				}
			}
		}

		for (int i = 0; i < line_length; ++i) {
			double h_k_0 = matrix_k_0[i][i];
			double h_k_2 = matrix_k_2[i][i];
			for (int j = i; j < line_length + 1; ++j) {
				matrix_k_0[i][j] /= h_k_0;
				matrix_k_2[i][j] /= h_k_2;
			}
		}
		std::cout << std::endl << std::endl;


		// back Gauss method
		double x_k_0 = 0, x_k_2 = 0;
		std::vector<double> vals_k_0(line_length), vals_k_2(line_length);
		for (int i = line_length - 1; i >= 0; --i) {
			if (i > line_length) {
				break;
			}

			x_k_0 = matrix_k_0[i][line_length];
			x_k_2 = matrix_k_2[i][line_length];
			for (int j = line_length - 1; j > i; --j) {
				if (i == line_length - 1 && j == line_length - 1) {
					vals_k_0[i] = matrix_k_0[i][j];
					vals_k_2[i] = matrix_k_2[i][j];
					continue;
				}
				else {
					x_k_0 -= matrix_k_0[i][j] * vals_k_0[j];
					x_k_2 -= matrix_k_2[i][j] * vals_k_2[j];
				}
			}
			vals_k_0[i] = x_k_0;
			vals_k_2[i] = x_k_2;
		}

		return { vals_k_0, vals_k_2 };
	}

	std::vector<double> GetCorrectAnswers() { return x_; }

	std::vector<double> GetVectorBK0() { return b_k_0_; }

	std::vector<double> GetVectorBK2() { return b_k_2_; }

	std::vector<std::vector<double>> GetMatrixK0() { return matrix_k_0_; }

	std::vector<std::vector<double>> GetMatrixK2() { return matrix_k_2_; }

private:
	std::vector<double> x_, b_k_0_, b_k_2_;
	std::vector<std::vector<double>> matrix_k_2_, matrix_k_0_;

};

int main() {
	Matrix m;
	std::cout << "Matrix with k = 0:" << std::endl;
	PrintMatrix(m.GetMatrixK0());
	std::cout << std::endl << std::endl << "Matrix with k = 2:" << std::endl;
	PrintMatrix(m.GetMatrixK2());

	auto answers = m.GaussWithoutMainElement();
	std::cout << std::endl << std::endl << "Answers for Gauss method without selection of the main element" << std::endl;
	std::cout << "Vector of answers for k = 0: " << std::endl;
	for (const auto& x : answers.first) {
		std::cout << x << " ";
	}
	std::cout << std::endl << std::endl << "Vector of answers for k = 2: " << std::endl;
	for (const auto& x : answers.second) {
		std::cout << x << " ";
	}
	std::cout << std::endl << std::endl;
	std::cout << "Error k = 0:" << CountInfelicity(m.GetCorrectAnswers(), answers.first) << std::endl;
	std::cout << "Error k = 2:" << CountInfelicity(m.GetCorrectAnswers(), answers.first) << std::endl;

	auto answers_with_main = m.GaussWithMainElement();
	std::cout << std::endl << "Answers for Gauss method with selection of the main element" << std::endl;
	std::cout << "Vector of answers for k = 0: " << std::endl;
	for (const auto& x : answers_with_main.first) {
		std::cout << x << " ";
	}
	std::cout << std::endl << std::endl << "Vector of answers for k = 2: " << std::endl;
	for (const auto& x : answers_with_main.second) {
		std::cout << x << " ";
	}
	std::cout << std::endl << std::endl;
	std::cout << "Error k = 0:" << CountInfelicity(m.GetCorrectAnswers(), answers_with_main.first) << std::endl;
	std::cout << "Error k = 2:" << CountInfelicity(m.GetCorrectAnswers(), answers_with_main.first) << std::endl;

	return 0;
}