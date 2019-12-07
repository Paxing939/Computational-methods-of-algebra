#include <cmath>
#include <iomanip>
#include "Functions.h"

double RoundNothing(double x) { return x; }

class Matrix {
public:

	~Matrix() {
		matrix_.clear();
	}

	// создание матрицы А
	Matrix() {
		srand(time(NULL));
		// случайное количество элементов матрицы
		int n = rand() % 3 + 10;
		// запоминаем числа для не диагональных элементов от -4 до 0
		std::vector<int> no_diag;
		for (int i = 0; i < (n * n - n) / 2; i++) {
			no_diag.push_back(rand() % 5 - 4);
		}

		matrix_.resize(n);
		for (auto& line : matrix_) {
			line.resize(n, -100);
		}

		// заполняем матрицу А
		int p = 0, sum_1_line = 0;
		double sum_1_col = 0.0f;
		for (int i = 0; i < matrix_.size(); i++) {
			for (int j = 0; j < matrix_.size(); j++) {
				if (i != j) {
					if (matrix_[i][j] == -100) {
						matrix_[i][j] = no_diag[p];
						matrix_[j][i] = no_diag[p];
						p++;
					}
				}
				if (i == 0 && j != 0) {
					sum_1_line += matrix_[0][j];
				}
				if (j == 0 && i != 0) {
					sum_1_col += matrix_[i][0];
				}
			}
		}

		matrix_k_0_ = matrix_;

		double sum_1_col_k_0 = (double)sum_1_line;
		matrix_k_0_[0][0] = -(double)sum_1_line + 1;
		sum_1_col_k_0 += matrix_k_0_[0][0];
		
		matrix_[0][0] = -(double)sum_1_line + pow(10, -2);
		sum_1_col += matrix_[0][0];

		for (int i = 1; i < matrix_.size(); i++) {
			matrix_k_0_[i][i] = sum_1_col_k_0 - matrix_k_0_[i][0];
			matrix_[i][i] = sum_1_col - matrix_[i][0];
		}

		int m = 17;
		for (int i = 0; i < n; i++) {
			x_.push_back(m + i);
		}
		// multiply x_ on matrix to get vector b
		// if after Gauss method answer vector and x_ will be the same
		// then we will know method works correctly
		b_ = MultMat(matrix_, x_);
		b_k_0_ = MultMat(matrix_k_0_, x_);

		// write matrices in the logs.txt
		LogMatrix(matrix_);
		LogMatrix(matrix_k_0_);
	}

	std::pair<std::vector<double>, std::vector<double>> Solve() {
		int line_length = matrix_.size();

		std::vector<std::vector<double>> matrix_transform(line_length);
		std::vector<std::vector<double>> matrix_transform_k_0(line_length);
		for (int i = 0; i < line_length; ++i) {
			matrix_transform[i].resize(line_length + 1);
			matrix_transform_k_0[i].resize(line_length + 1);
			for (int j = 0; j < line_length + 1; j++) {
				if (j != line_length) {
					matrix_transform[i][j] = matrix_[i][j];
					matrix_transform_k_0[i][j] = matrix_k_0_[i][j];
				}
				else {
					matrix_transform[i][j] = b_[i];
					matrix_transform_k_0[i][j] = b_k_0_[i];
				}
			}
		}

		// forward
		D_.resize(line_length);
		D_k_0_.resize(line_length);
		for (int k = 0; k < line_length; k++) {
			D_[k].resize(line_length);
			D_[k][k] = matrix_transform[k][k];
			D_k_0_[k].resize(line_length);
			D_k_0_[k][k] = matrix_transform_k_0[k][k];
			std::vector<double> t(line_length), t0(line_length);
			for (int i = k + 1; i < line_length; i++) {
				t0[i] = matrix_transform_k_0[i][k];
				matrix_transform_k_0[i][k] = matrix_transform_k_0[i][k] / matrix_transform_k_0[k][k];

				t[i] = matrix_transform[i][k];
				matrix_transform[i][k] = matrix_transform[i][k] / matrix_transform[k][k];

				for (int j = k + 1; j <= i; j++) {
					matrix_transform[i][j] -= matrix_transform[i][k] * t[j];
					matrix_transform_k_0[i][j] -= matrix_transform_k_0[i][k] * t0[j];
				}
			}
		}

		// L and D
		std::vector<std::vector<double>> Lt(matrix_.size()), Lt_k_0(matrix_k_0_.size());
		L_.resize(matrix_.size());
		L_k_0_.resize(matrix_.size());
		for (int i = 0; i < line_length; ++i) {
			for (int j = 0; j < line_length; ++j) {
				if (i > j) {
					L_[i].push_back(matrix_transform[i][j]);
					L_k_0_[i].push_back(matrix_transform_k_0[i][j]);
				} else {
					if (i == j) {
						L_[i].push_back(1.0f);
						L_k_0_[i].push_back(1.0f);
					} else {
						L_[i].push_back(0.0f);
						L_k_0_[i].push_back(0.0f);
					}
				}
			}
		}

		// transpose L
		for (int i = 0; i < line_length; ++i) {
			Lt[i].resize(line_length, 0);
			Lt_k_0[i].resize(line_length, 0);
			for (int j = 0; j < line_length; ++j) {
				if (i < j) {
					Lt[i][j] = L_[j][i];
					Lt_k_0[i][j] = L_k_0_[j][i];
				} else {
					if (i == j) {
						Lt[i][j] = 1;
						Lt_k_0[i][j] = 1;
					}
				}
			}
		}
 
		std::cout << "Matrixs L_ and Lt with k = 2:";
		PrintMatrix(L_); PrintMatrix(Lt);
		std::cout << std::endl << "Matrixs L_ and Lt with k = 0:";
		PrintMatrix(L_k_0_); PrintMatrix(Lt_k_0);

		// back and solve Ly = b
		std::vector<double> y(L_.size()), y_k_0(L_k_0_.size());
		for (int i = 0; i < L_.size(); i++) {
			y[i] = b_[i];
			y_k_0[i] = b_k_0_[i];
			for (int j = i - 1; j >= 0; j--) {
				y[i] -= L_[i][j] * y[j];
				y_k_0[i] -= L_k_0_[i][j] * y_k_0[j];
			}
		}

		// solve Dz = y
		for (int i = 0; i < L_.size(); i++) {
			y[i] /= D_[i][i];
			y_k_0[i] /= D_k_0_[i][i];
		}

		// solve Ltx = z
		std::vector<double> answer(y.size()), answer_k_0(y_k_0.size());
		for (int i = y.size() - 1; i >= 0; i--) {
			answer[i] = y[i];
			answer_k_0[i] = y_k_0[i];
			for (int j = i + 1; j < L_.size(); j++) {
				answer[i] -= Lt[i][j] * answer[j];
				answer_k_0[i] -= Lt_k_0[i][j] * answer_k_0[j];
			}
		}

		LogVector(answer);
		LogVector(answer_k_0);
		LogMatrix(L_);
		LogMatrix(L_k_0_);
		std::vector<double> d, d_k_0;
		for (int i = 0; i < D_.size(); i++) {
			d.push_back(D_[i][i]);
			d_k_0.push_back(D_k_0_[i][i]);
		}
		LogVector(d);
		LogVector(d_k_0);
		for (int i = 0; i < line_length; ++i) {
			for (int j = 0; j < line_length; ++j) {
				std::cout << std::setw(7) << std::setprecision(2) << std::left << matrix_transform[i][j] << " ";
			}
			std::cout << " " << matrix_transform[i][line_length] << std::endl;
		}

		std::cout << "Errors:" << std::endl << "k = 2:" << CountInfelicity(x_, answer) << std::endl << "k = 0:" << CountInfelicity(x_, answer_k_0);
		std::ofstream fout("logs.txt", std::ios::app);
		fout << "k = 2:" << CountInfelicity(x_, answer) << std::endl << "k = 0:" << CountInfelicity(x_, answer_k_0);
		fout.close();
		return std::pair<std::vector<double>, std::vector<double>> {answer, answer_k_0};
	}

private:
	std::vector<std::vector<double>> matrix_, L_, D_;
	std::vector<std::vector<double>> matrix_k_0_, L_k_0_, D_k_0_;
	std::vector<double> b_, x_, b_k_0_;
};

int main() {
	remove("logs.txt");
	Matrix m;
	std::cout << std::endl;
	std::pair<std::vector<double>, std::vector<double>> pair = m.Solve();
	std::cout << std::endl << "Vector with k = -2" << std::endl;
	for (const auto& el : pair.first) {
		std::cout << std::setw(7) << std::setprecision(7) << el << " ";
	}
	std::cout << std::endl << "Vector with k = 0" << std::endl;
	for (const auto& el : pair.second) {
		std::cout << std::setw(7) << std::setprecision(5) << el << " ";
	}
	return 0;
}