#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <fstream>

//double R(double x) {
//	int k = 1000;
//	return round(x * k) / (double)k;
//}
double R(double x) { return x; }

double CountInfelicity(std::vector<double> x_de_ure, std::vector<double> x_de_facto) {
	double sum_x2 = 0, sum_x = 0;
	for (int i = 0; i < x_de_facto.size(); i++) {
		sum_x2 += x_de_facto[i] * x_de_facto[i];
		sum_x += (x_de_facto[i] - x_de_ure[i]) * (x_de_facto[i] - x_de_ure[i]);
	}
	return sum_x / sum_x2;
}

void LogMatrix(std::vector<std::vector<double>> matrix) {
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

void LogVector(std::vector<double> vector) {
	std::ofstream fout("logs.txt", std::ios::app);
	for (const auto& el : vector) {
		fout << el << '\t';
	}
	fout << std::endl << std::endl;
	fout.close();
}

std::vector<std::vector<double>> MultMatrix(std::vector<std::vector<double>> a, std::vector<std::vector<double>> b_) {
	std::vector<std::vector<double>> c(a.size());
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

void PrintMatrix(const std::vector<std::vector<double>>& matrix) {
	for (const std::vector<double>& vec : matrix) {
		for (const double& x_ : vec) {
			std::cout << std::setw(7) << std::setprecision(2) << std::left << x_ << " ";
		}
		std::cout << std::endl;
	}
}

class Matrix {
public:

	~Matrix() {
		matrix_.clear();
	}

	// создание матрицы ј
	Matrix() {
		srand(time(NULL));
		// случайное количество элементов матрицы
		int n = rand() % 3 + 10;
		// запоминаем числа дл€ не диагональных элементов от -4 до 0
		std::vector<int> no_diag;
		for (int i = 0; i < (n * n - n) / 2; i++) {
			no_diag.push_back(rand() % 5 - 4);
		}

		matrix_.resize(n);
		for (auto& line : matrix_) {
			line.resize(n, -100);
		}

		// заполн€ем матрицу ј
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

		// переписываем матрицу ј с к = 0
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
		// с помощью умножени€ матриц находим вектора b
		b_ = MultMat(x_, matrix_);
		b_k_0_ = MultMat(x_, matrix_k_0_);
		// пишем исходные матрицы в файл
		LogMatrix(matrix_);
		LogMatrix(matrix_k_0_);
	}

	// функци€ умножени€ матриц
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

	// функци€ решени€
	std::pair<std::vector<double>, std::vector<double>> Solve() {
		int line_length = matrix_.size();
		// запол€нем матрицу гаусса
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

		// пр€мой ход гаусса
		D_.resize(line_length);
		D_k_0_.resize(line_length);
		for (int k = 0; k < line_length; k++) {
			D_[k].resize(line_length);
			D_[k][k] = matrix_transform[k][k];
			D_k_0_[k].resize(line_length);
			D_k_0_[k][k] = matrix_transform_k_0[k][k];

			for (int i = k + 1; i < line_length; i++) {
				double t_k_0 = R(matrix_transform_k_0[i][k] / matrix_transform_k_0[k][k]);
				matrix_transform_k_0[i][k] = t_k_0;
				double t = R(matrix_transform[i][k] / matrix_transform[k][k]);
				matrix_transform[i][k] = t;

				for (int j = k + 1; j < line_length; ++j) {
					matrix_transform[i][j] -= R(matrix_transform[k][j] * t);
					matrix_transform_k_0[i][j] -= R(matrix_transform_k_0[k][j] * t_k_0);
				}
			}
		}

		// заполн€ем L из D
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

		// транспонируем матрицу L
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

		// выводим в консоль L 
		std::cout << "Matrixs L_ and Lt with k = 2:";
		PrintMatrix(L_); PrintMatrix(Lt);
		std::cout << std::endl << "Matrixs L_ and Lt with k = 0:";
		PrintMatrix(L_k_0_); PrintMatrix(Lt_k_0);

		// обратный ход гаусса
		// решаем уравнение Ly = b
		std::vector<double> y(L_.size()), y_k_0(L_k_0_.size());
		for (int i = 0; i < L_.size(); i++) {
			y[i] = b_[i];
			y_k_0[i] = b_k_0_[i];
			for (int j = i - 1; j >= 0; j--) {
				y[i] -= R(L_[i][j] * y[j]);
				y_k_0[i] -= R(L_k_0_[i][j] * y_k_0[j]);
			}
		}

		// решаем уравнение Dz = y
		for (int i = 0; i < L_.size(); i++) {
			y[i] /= R(D_[i][i]);
			y_k_0[i] /= R(D_k_0_[i][i]);
		}

		// решаем уравнение Ltx = z
		std::vector<double> answer(y.size()), answer_k_0(y_k_0.size());
		for (int i = y.size() - 1; i >= 0; i--) {
			answer[i] = y[i];
			answer_k_0[i] = y_k_0[i];
			for (int j = i + 1; j < L_.size(); j++) {
				answer[i] -= R(Lt[i][j] * answer[j]);
				answer_k_0[i] -= R(Lt_k_0[i][j] * answer_k_0[j]);
			}
		}

		// пишем в файл ответы
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
		// считаем погрешности
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