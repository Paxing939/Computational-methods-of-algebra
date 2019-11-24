#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <fstream>

float R(float x) {
	int k = 1000;
	return round(x * k) / (float)k;
}

float CountInfelicity(std::vector<float> x_de_ure, std::vector<float> x_de_facto) {
	float sum_x2 = 0, sum_x = 0;
	for (int i = 0; i < x_de_facto.size(); i++) {
		sum_x2 += x_de_facto[i] * x_de_facto[i];
		sum_x += (x_de_facto[i] - x_de_ure[i]) * (x_de_facto[i] - x_de_ure[i]);
	}
	std::ofstream fout("errors.txt");
	fout << sum_x / sum_x2 << std::endl;
	fout.close();
	return sum_x / sum_x2;
}

void LogMatrix(std::vector<std::vector<float>> matrix) {
	std::ofstream fout("logs.txt", std::ios::app);
	for (const auto& vec : matrix) {
		for (auto& el : vec) {
			fout << el << '\t';
		}
		fout << std::endl;
	}
	fout << std::endl << std::endl;
	fout.close();
}

void LogVector(std::vector<float> vector) {
	std::ofstream fout("logs.txt", std::ios::app);
	for (const auto& el : vector) {
		fout << el << '\t';
	}
	fout << std::endl << std::endl;
	fout.close();
}

std::vector<std::vector<float>> MultMatrix(std::vector<std::vector<float>> a, std::vector<std::vector<float>> b_) {
	std::vector<std::vector<float>> c(a.size());
	for (int i = 0; i < a.size(); i++) {
		c[i].resize(b_[i].size());
		for (int j = 0; j < b_[i].size(); j++) {
			c[i][j] = 0;
			for (int k = 0; k < a[i].size(); k++)
				c[i][j] += a[i][k] * b_[k][j];
		}
	}
	return c;
}

void PrintMatrix(const std::vector<std::vector<float>>& matrix) {
	for (const auto& vec : matrix) {
		for (const auto& x_ : vec) {
			std::cout << std::setw(7) << std::setprecision(2) << std::left << x_ << " ";
		}
		std::cout << std::endl;
	}
}

void PrintVector(const std::vector<float>& vector) {
	for (const auto& x_ : vector) {
		std::cout << std::setw(7) << std::setprecision(2) << std::left << x_ << " ";
	}
	std::cout << std::endl;
}

class Matrix {
public:

	// функция умножения матриц
	std::vector<float> MultMat(std::vector<float> a, std::vector<std::vector<float>> matrix) {
		std::vector<float> c(matrix.size());

		for (int i = 0; i < matrix.size(); i++) {
			c[i] = 0;
			for (int k = 0; k < matrix.size(); k++) {
				c[i] += matrix[i][k] * a[k];
			}
		}
		return c;
	}

	~Matrix() {
		matrix_.clear();
	}

	// создание матрицы А
	Matrix() {
		srand(time(NULL));
		// случайное количество элементов матрицы
		int n = 4;// rand() % 3 + 10;

		matrix_.resize(n, std::vector<float>(n));
		for (auto& line : matrix_) {
			line.resize(n);
		}

		// заполняем матрицу А
		// заполнение недиагональных элементов
		int p = 0, sum_1_col = 0;
		for (int i = 0; i < matrix_.size(); i++) {
			for (int j = 0; j < matrix_.size(); j++) {
				matrix_[i][j] = rand() % 100 - 50;	
			}
		}
		// пишем исходную матрицу в файл
		LogMatrix(matrix_);
		PrintMatrix(matrix_);
		std::cout << std::endl;
	}

	void ImplementDanilevski() {
		std::vector<std::vector<float>> M_i, M_i_1;
		while (true) {
			for (int i = matrix_.size() - 2; i >= 0; i--) {
				std::cout << i << std::endl;
				if (i == 0) {
					int sfgvds = 1010234;
				}
				M_i = CreateM_i(i);
				std::cout << std::endl;
				PrintMatrix(M_i);
				LogMatrix(M_i);

				std::ofstream fout("logs.txt", std::ios::app);
				fout << std::endl << std::endl;
				fout.close();

				M_i_1 = CreateM_i_1(i);
				std::cout << std::endl;
				PrintMatrix(M_i_1);
				LogMatrix(M_i_1);

				fout.open("logs.txt", std::ios::app);
				fout << std::endl << std::endl;
				fout.close();

				matrix_ = MultMatrix(matrix_, M_i);
				matrix_ = MultMatrix(M_i_1, matrix_);
				PrintMatrix(matrix_);
				LogMatrix(matrix_);


				fout.open("logs.txt", std::ios::app);
				fout << std::endl << std::endl << "----------------------------------";
				fout.close();

				if (i != 0) {
					if (matrix_[i][i - 1] == 0) {
						continue;
					}
				}
			}
			break;
		}
	}

	std::vector<std::vector<float>> CreateM_i(int n) {
		std::vector<std::vector<float>> M_i(matrix_.size(), std::vector<float>(matrix_.size()));

		for (int i = 0; i < matrix_.size(); i++) {
			if (i != n) {
				M_i[i][i] = 1;
			}
			if (i == n) {
				M_i[n][i] = 1 / matrix_[n + 1][n];
			} else {
				M_i[n][i] = -matrix_[n + 1][i] / matrix_[n + 1][n];
			}
		}

		return M_i;
	}

	std::vector<std::vector<float>> CreateM_i_1(int n) {
		std::vector<std::vector<float>> M_i_1(matrix_.size(), std::vector<float>(matrix_.size()));

		for (int i = 0; i < matrix_.size(); i++) {
			if (i != n) {
				M_i_1[i][i] = 1;
			}
			M_i_1[n][i] = matrix_[n + 1][i];
		}

		return M_i_1;
	}

private:
	std::vector<std::vector<float>> matrix_;
	std::vector<float> f_, x_, x_curr_;
	int k_max_, k_;
	float E_;
};

int main() {
	remove("logs.txt");
	Matrix m;
	std::cout << std::endl;
	m.ImplementDanilevski();

	return 0;
}