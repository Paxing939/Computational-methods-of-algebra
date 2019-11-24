#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <fstream>

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

void PrintMatrix(const std::vector<std::vector<float>>& matrix) {
	for (const std::vector<float>& vec : matrix) {
		for (const float& x_ : vec) {
			std::cout << std::setw(7) << std::setprecision(2) << std::left << x_ << " ";
		}
		std::cout << std::endl;
	}
}

class Matrix {

public:
	Matrix() { // создание матрицы для 1 случая
	//	line_length = (15 + rand() % 5);
	//	matrix_.clear();
	//	matrix_.resize(line_length);
	//	for (auto& vec : matrix_) {
	//		vec.resize(line_length);
	//	}
	//	for (auto& vec : matrix_) {
	//		for (auto& x : vec) {
	//			x = rand() % 200 - 100;
	//		}
	//	}


		srand(time(NULL));
		// случайное количество элементов матрицы
		int n = 5;//rand() % 3 + 12;
		// запоминаем числа для не диагональных элементов от -4 до 0
		std::vector<int> no_diag;
		for (int i = 0; i < (n * n - n) / 2; i++) {
			no_diag.push_back(rand() % 5 - 4);
		}

		matrix_k_2_.resize(n);
		for (auto& line : matrix_k_2_) {
			line.resize(n, -100);
		}

		// заполняем матрицу А
		int p = 0;
		float sum_1_col = 0.0f, sum_1_line = 0.0f;
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

		// переписываем матрицу А с к = 0
		matrix_k_0_ = matrix_k_2_;

		matrix_k_0_[0][0] = -sum_1_line + 1;
		matrix_k_2_[0][0] = -sum_1_line + pow(10, -2);

		float sum_1_col_k_0 = 0.0f, sum_1_col_k_2 = 0.0f;
		for (int i = 1; i < matrix_k_0_.size(); i++) {
			sum_1_col_k_0 = 0.0f;
			sum_1_col_k_2 = 0.0f;
			for (int j = 0; j < line_length; j++) {
				if (i != j) {
					sum_1_col_k_0 += matrix_k_0_[i][j];
					sum_1_col_k_2 += matrix_k_2_[i][j];
				}
			}
			matrix_k_0_[i][i] = sum_1_col - matrix_k_0_[i][0];
		}
		PrintMatrix(matrix_k_0_);
		std::cout << std::endl;
		PrintMatrix(matrix_k_2_);
		//int m = 17;
		//for (int i = 0; i < n; i++) {
		//	x_.push_back(m + i);
		//}
		//// с помощью умножения матриц находим вектора b
		//b_ = MultMat(x_, matrix_);
		//b_k_0_ = MultMat(x_, matrix_k_0_);
		//// пишем исходные матрицы в файл
		//LogMatrix(matrix_);
		//LogMatrix(matrix_k_0_);
	}


	Matrix(const int& n1) { // создание матрицы для 2 случая
	//	line_length = n1;
	//	matrix_.resize(line_length);
	//	for (int i = 0; i < line_length; ++i) {
	//		matrix_[i].resize(line_length);
	//	}
	//	for (int i = 0; i < line_length; ++i) {
	//		for (int j = 0; j < line_length; ++j) {
	//			if (i != j) {
	//				matrix_[i][j] = (-1) * (rand() % 4);
	//			}
	//			else {
	//				matrix_[i][j] = 0;
	//				if (i == 0) {
	//					matrix_[i][j] = pow(10, -2);
	//				}
	//			}
	//		}
	//	}
	//	for (int i = 0; i < line_length; ++i) {
	//		for (int j = 0; j < line_length; ++j) {
	//			if (i != j) {
	//				matrix_[i][i] += -matrix_[i][j];
	//			}
	//		}
	//	}
	//	PrintMatr();
	//	std::cout << std::endl << std::endl;
	//}

	//~Matrix() {
	//	matrix_.clear();
	}

	std::vector<float> GaussWithoutMainElement() {
		// заполнение матрицы
	//	std::vector<float> vec(line_length);
	//	int m = 2;
	//	for (int i = 0; i < line_length; ++i) {
	//		vec[i] = m;
	//		++m;
	//	}
	//	std::vector<float> vec1(line_length);
	//	for (int i = 0; i < line_length; ++i) {
	//		for (int j = 0; j < line_length; ++j) {
	//			vec1[i] += matrix_[i][j] * vec[j];
	//		}
	//	}
	//	// заполняем матрицу гаусса
	//	std::vector<std::vector<float>> matrix_transform(line_length);
	//	for (int i = 0; i < line_length; ++i) {
	//		matrix_transform[i].resize(line_length + 1);
	//	}
	//	for (int i = 0; i < line_length; ++i) {
	//		for (int j = 0; j < line_length + 1; ++j) {
	//			if (j != line_length) {
	//				matrix_transform[i][j] = matrix_[i][j];
	//			}
	//			else {
	//				matrix_transform[i][j] = vec1[i];
	//			}
	//		}
	//	}

	//	// прямой ход гаусса
	//	for (int i = 0; i < line_length; ++i) {
	//		for (int s = i + 1; s < line_length; ++s) {
	//			float h = matrix_transform[s][i] / matrix_transform[i][i];
	//			for (int j = i; j < line_length + 1; ++j) {
	//				matrix_transform[s][j] -= matrix_transform[i][j] * h;
	//			}
	//		}
	//	}

	//	// выводим на экран
	//	for (int i = 0; i < line_length; ++i) {
	//		for (int j = 0; j < line_length; ++j) {
	//			std::cout << std::setw(7) << std::setprecision(2) << std::left << matrix_transform[i][j] << " ";
	//		}
	//		std::cout << " " << matrix_transform[i][line_length] << std::endl;
	//	}

	//	// вычисляем ответ
	//	std::vector<float> vals(line_length);
	//	float x = 0;
	//	for (int i = line_length - 1; i >= 0; --i) {
	//		if (i > line_length) {
	//			break;
	//		}
	//		x = matrix_transform[i][line_length];
	//		for (int j = line_length - 1; j > i; --j) {
	//			if (i == line_length - 1 && j == line_length - 1) {
	//				vals[i] = matrix_transform[i][j];
	//				continue;
	//			}
	//			else {
	//				x -= matrix_transform[i][j] * vals[j];
	//			}
	//		}
	//		vals[i] = x / matrix_transform[i][i];
	//	}

	//	std::cout << "Error :" << CountInfelicity(vec, vals) << std::endl;
	//	return vals;
		return std::vector<float>{};
	}

	std::vector<float> GaussWithMainElement() {
	//	std::vector<float> vals(line_length), vec(line_length), vec1(line_length);
	//	int m = 2;
	//	for (int i = 0; i < line_length; ++i) {
	//		vec[i] = m;
	//		++m;
	//	}
	//	for (int i = 0; i < line_length; ++i) {
	//		for (int j = 0; j < line_length; ++j) {
	//			vec1[i] += matrix_[i][j] * vec[j];
	//		}
	//	}

	//	std::vector<std::vector<float>> matrix_transform(line_length);
	//	for (int i = 0; i < line_length; ++i) {
	//		matrix_transform[i].resize(line_length + 1);
	//	}
	//	for (int i = 0; i < line_length; ++i) {
	//		for (int j = 0; j < line_length + 1; ++j) {
	//			if (j != line_length) {
	//				matrix_transform[i][j] = matrix_[i][j];
	//			} else {
	//				matrix_transform[i][j] = vec1[i];
	//			}
	//		}
	//	}

	//	for (int i = 0; i < line_length; ++i) {
	//		for (int j = 0; j < line_length; ++j) {
	//			std::cout << std::setw(6) << std::setprecision(4) << std::left << matrix_transform[i][j] << " ";
	//		}
	//		std::cout << " " << matrix_transform[i][line_length] << std::endl;
	//	}

	//	int max_index = 0;
	//	float max = 0;
	//	for (int k = 0; k < line_length; ++k) {
	//		max_index = k;
	//		max = 0;
	//		for (int i = k; i < line_length; ++i) {
	//			if (matrix_transform[i][k] > max) {
	//				max = matrix_transform[i][k];
	//				max_index = i;
	//			}
	//		}
	//		swap(matrix_transform[k], matrix_transform[max_index]);
	//	}
	//	for (int i = 0; i < line_length; ++i) {
	//		for (int s = i + 1; s < line_length; ++s) {
	//			float h = matrix_transform[s][i] / matrix_transform[i][i];
	//			for (int j = i; j < line_length + 1; ++j) {
	//				matrix_transform[s][j] -= matrix_transform[i][j] * h;
	//			}
	//		}
	//	}
	//	for (int i = 0; i < line_length; ++i) {
	//		float h = matrix_transform[i][i];
	//		for (int j = i; j < line_length + 1; ++j) {
	//			matrix_transform[i][j] /= h;
	//		}
	//	}
	//	std::cout << std::endl << std::endl;

	//	for (int i = 0; i < line_length; ++i) {
	//		for (int j = 0; j < line_length; ++j) {
	//			std::cout << std::setw(6) << std::setprecision(4) << std::left << matrix_transform[i][j] << " ";
	//		}
	//		std::cout << " " << matrix_transform[i][line_length] << std::endl;
	//	}
	//	float x = 0;
	//	for (int i = line_length - 1; i >= 0; --i) {
	//		if (i > line_length) {
	//			break;
	//		}
	//		x = matrix_transform[i][line_length];
	//		for (int j = line_length - 1; j > i; --j) {
	//			if (i == line_length - 1 && j == line_length - 1) {
	//				vals[i] = matrix_transform[i][j];
	//				continue;
	//			}
	//			else {
	//				x -= matrix_transform[i][j] * vals[j];
	//			}
	//		}
	//		vals[i] = x;
	//	}
	//	std::cout << "Error :" << CountInfelicity(vec, vals) << std::endl;
		return std::vector<float>{};//vals;
	}

	//void PrintMatr() {
	//	for (const auto& vec : matrix_) {
	//		for (const auto& x : vec) {
	//			std::cout << std::left << std::setw(7) << x << " ";
	//		}
	//		std::cout << std::endl;
	//	}
	//}

private:
	int line_length;
	std::vector<std::vector<float>> matrix_k_2_, matrix_k_0_;

};

int main() {
	Matrix m;
	std::vector<float> vec = m.GaussWithoutMainElement();
	std::cout << std::endl << std::endl;
	for (const auto& x : vec) {
		std::cout << x << " ";
	}
	std::cout << std::endl << std::endl;
	//Matrix x;
	std::vector<float> vec2 = m.GaussWithMainElement();
	for (const auto& x : vec2) {
		std::cout << x << " ";
	}
	return 0;
}