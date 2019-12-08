#include <cmath>
#include <iomanip>
#include "Functions.h"

float Norm(const std::vector<float>& vec) {
	float max = 0;
	for (const auto& el : vec) {
		if (max < abs(el)) {
			max = abs(el);
		}
	}
	return max;
}

void Normalization(std::vector<float>& vec) {
	float max = 0;
	for (const auto& el : vec) {
		if (el > max) {
			max = el;
		}
	}
	for (auto& el : vec) {
		el /= max;
	}
}

float ScalarProduct(const std::vector<float>& vec_1, const std::vector<float>& vec_2) {
	float ret_val = 0;
	for (int i = 0; i < vec_1.size(); i++) {
		ret_val += vec_1[i] * vec_2[i];
	}
	return ret_val;
}

float Sign(float number) {
	if (number < 0) {
		return -1;
	}
	if (number > 0) {
		return 1;
	}
	return 0;
}

std::vector<float> Eigenvectors(const std::vector<float>& u, const std::vector<float>& v, float lambda) {
	std::vector<float> eigenvector(u.size());
	for (int i = 0; i < v.size(); i++) {
		eigenvector[i] = v[i] - lambda * u[i];
	}
	return eigenvector;
}

float CountSecondEigenValue(std::vector<float> v_m, std::vector<float> v_m1, 
								std::vector<float> u_m, std::vector<float> u_m1, float lambda) {
	float max = 0;
	int index = 0;
	for(int i = 0; i < v_m.size(); i++) {
		if (abs(v_m[i] - lambda * u_m1[i]) > max) {
			max = abs(v_m[i] - lambda * u_m1[i]);
			index= i;
		}
	}
	return (v_m1[index] * Norm(v_m) - lambda * v_m[index]) / (v_m[index] - lambda * u_m1[index]);
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
		int n = 4;//rand() % 3 + 10;
		// запоминаем числа дл€ не диагональных элементов от -4 до 0
		std::vector<int> no_diag;
		for (int i = 0; i < (n * n - n) / 2; i++) {
			no_diag.push_back(rand() % 5 - 4);
		}

		matrix_.resize(n);
		for (auto& line : matrix_) {
			line.resize(n, -6);
		}

		// fill matrix_
		int p = 0, sum_1_line = 0;
		float sum_1_col = 0.0f;
		for (int i = 0; i < matrix_.size(); i++) {
			for (int j = 0; j < matrix_.size(); j++) {
				if (i != j) {
					if (matrix_[i][j] == -6) {
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

		matrix_[0][0] = -(float)sum_1_line + pow(10, -2);
		sum_1_col += matrix_[0][0];

		for (int i = 1; i < matrix_.size(); i++) {
			matrix_[i][i] = sum_1_col - matrix_[i][0];
		}
		matrix_ = { { -23, -18, -35, 28 },
					{ 5, 10, 5, -9 },
					{ -49, 12, 25, 16 },
					{ 26, -48, -42, 21} };
		// write matrices in the logs.txt
		LogMatrix(matrix_);
	}

	void StartIterations() {
		int line_length = matrix_.size();
		std::vector<float> x_1(line_length), x_2(line_length),
			v_30, v_31, u_30, u_29, v_49, v_50, u_49, u_48;
		x_2[0] = 1;
		float max = 0;
		for (int i = 0; i < 46; i++) {
			if (i == 30) {
				v_30 = x_2;
			}
			if (i == 31) {
				v_31 = x_2;
			}
			Normalization(x_2);
			if (i == 29) {
				u_29 = x_2;
			}
			if (i == 30) {
				u_30 = x_2;
			}
			std::cout << "k = " << i << ":   ";
			PrintVector(x_2);
			x_2 = MultMat(matrix_, x_2);
			max = 0;
		}
		std::cout << std::endl;
		float lamda_simple = 0, lamda_scalar_product = 0;
		for (int i = 46; i < 51; i++) {
			if (i == 49) {
				v_49 = x_2;
			}
			if (i == 50) {
				v_50 = x_2;
			}
			x_1 = x_2;
			Normalization(x_1);
			if (i == 49) {
				u_49 = x_1;
			}
			if (i == 48) {
				u_48 = x_1;
			}
			std::cout << "k = " << i << ":   ";
			PrintVector(x_1);
			x_2 = MultMat(matrix_, x_1);
			lamda_simple = Norm(x_2) * Sign(Norm(x_1));
			lamda_scalar_product = ScalarProduct(x_2, x_1) / ScalarProduct(x_1, x_1);
		}
		std::vector<float> eigenvector_simple = Eigenvectors(x_1, x_2, lamda_simple);
		std::vector<float> eigenvector_scalar = Eigenvectors(x_1, x_2, lamda_scalar_product);
		PrintVector(eigenvector_simple);
		std::cout << std::endl << std::endl;
		PrintVector(eigenvector_simple);
		std::cout << std::endl << std::endl;
		std::cout << "Eigenvalue by lambda = v / u: " << lamda_simple << std::endl
				<< "Eigenvalue by lambda <v, u> / <u, u>: " << lamda_scalar_product;

		float eigenvalue_30 = CountSecondEigenValue(v_30, v_31, u_30, u_29, lamda_simple);
		float eigenvalue_50 = CountSecondEigenValue(v_49, v_50, u_49, u_48, lamda_simple);
		float eigenvalue_50_scalar_eigenval = CountSecondEigenValue(v_49, v_50, u_49, u_48, lamda_scalar_product);
		std::cout << std::endl << std::endl << eigenvalue_30 << "  " << eigenvalue_50 << "  " << eigenvalue_50_scalar_eigenval << std::endl;


	}

private:
	std::vector<std::vector<float>> matrix_;
	std::vector<float> x_;
};

int main() {
	remove("logs.txt");
	Matrix m;
	std::cout << std::endl;
	m.StartIterations();

	return 0;
}