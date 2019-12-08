#include <cmath>
#include <iomanip>
#include "Functions.h"

class Matrix {
public:

	// функция умножения матриц
	std::vector<double> MultMat(std::vector<double> a, std::vector<std::vector<double>> matrix) {
		std::vector<double> c(matrix.size());

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
		int n = 4;//rand() % 3 + 310;

		matrix_.resize(n, std::vector<double>(n));
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

		std::ofstream fout("matrix.txt");
		fout << "{";
		for (int i = 0; i < matrix_.size(); i++) {
			fout << "{";
			for (int j = 0; j < matrix_.size() - 1; j++) {
				fout << matrix_[i][j] << ", ";
			}
			if (i == matrix_.size() - 1)
				fout << matrix_[i][matrix_.size() - 1] << "}";
			else
				fout << matrix_[i][matrix_.size() - 1] << "},";
		}
		fout << "}";
	}

	void ImplementDanilevski() {
		std::vector<std::vector<double>> M_i, M_i_1;
		while (true) {
			for (int i = matrix_.size() - 2; i >= 0; i--) {
				std::cout << "M_i" << std::endl;
				M_i = CreateM_i(i);
				std::cout << std::endl;
				PrintMatrix(M_i);
				LogMatrix(M_i);

				std::ofstream fout("logs.txt", std::ios::app);
				fout << std::endl << std::endl;
				fout.close();

				M_i_1 = CreateM_i_1(i);
				std::cout << "M_i_1" << std::endl;
				PrintMatrix(M_i_1);
				LogMatrix(M_i_1);

				fout.open("logs.txt", std::ios::app);
				fout << std::endl << std::endl;
				fout.close();

				matrix_ = MultMatrix(matrix_, M_i);
				std::cout << "M_i_1 multiply on A" << std::endl;
				PrintMatrix(matrix_);
				matrix_ = MultMatrix(M_i_1, matrix_);
				std::cout << std::endl << "----------------------------------" << std::endl;
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

	std::vector<std::vector<double>> CreateM_i(int n) {
		std::vector<std::vector<double>> M_i(matrix_.size(), std::vector<double>(matrix_.size()));

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

	std::vector<std::vector<double>> CreateM_i_1(int n) {
		std::vector<std::vector<double>> M_i_1(matrix_.size(), std::vector<double>(matrix_.size()));

		for (int i = 0; i < matrix_.size(); i++) {
			if (i != n) {
				M_i_1[i][i] = 1;
			}
			M_i_1[n][i] = matrix_[n + 1][i];
		}

		return M_i_1;
	}

private:
	std::vector<std::vector<double>> matrix_;
	std::vector<double> f_, x_, x_curr_;
	int k_max_, k_;
	double E_;
};

int main() {
	remove("logs.txt");
	Matrix m;
	std::cout << std::endl;
	m.ImplementDanilevski();

	return 0;
}