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

	~Matrix() {
		matrix_.clear();
	}

	// создание матрицы А
	Matrix() {
		srand(time(NULL));
		// случайное количество элементов матрицы
		int N = rand() % 3 + 10, m = 17, k = 2;
		matrix_.resize(N + 1);
		matrix_[0].resize(N + 1);
		matrix_[0][0] = m, matrix_[0][1] = m - 1;
		y_.push_back(1);

		// заполняем матрицу А
		for (int i = 1; i < N; i++)	{
			matrix_[i].resize(N + 1, 0.0f);
			matrix_[i][i - 1] = -k;
			matrix_[i][i] = m + k + i - 1;
			matrix_[i][i + 1] = m + i - 1;
			y_.push_back(i + 1);
		}
		matrix_[N].resize(N + 1, 0.0f);
		matrix_[N][N - 1] = -k;
		matrix_[N][N] = m + k + N - 1;
		y_.push_back(N + 1);

		f_ = MultMat(y_, matrix_);
	}

	std::vector<float> GetY() { return y_; }

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

	std::vector<float> Solve() {
		// вывод на экран и в файл
		LogMatrix(matrix_); LogVector(f_);
		std::cout << "Vector y we need to get: " << std::endl;
		PrintVector(y_);
		std::cout << "Input matrix: " << std::endl;
		PrintMatrix(matrix_);
		std::cout << std::endl;
		std::cout << "Vector f: " << std::endl;
		PrintVector(f_);
		std::cout << std::endl;
		std::cout << std::endl;

		// зануляем нижнюю диагональ
		matrix_[0][1] /= matrix_[0][0]; f_[0] /= matrix_[0][0];
		matrix_[0][0] = 1;
		for (int i = 1; i < matrix_.size(); i++) {
			float dev = R(matrix_[i][i - 1] / matrix_[i - 1][i - 1]);
			matrix_[i][i]		-= R(matrix_[i - 1][i] * dev);
			f_[i]				-= R(f_[i - 1] * dev);
			matrix_[i][i - 1]	-= R(matrix_[i - 1][i - 1] * dev);
		}

		// вывод на экран и в файл
		LogMatrix(matrix_); LogVector(f_);
		std::cout << "Matrix after first step: " << std::endl;
		PrintMatrix(matrix_);
		std::cout << "Vector f after first step: " << std::endl;
		PrintVector(f_);

		// вычисляем вектор решений
		std::vector<float> answers(matrix_.size());
		answers[answers.size() - 1] = R(f_[f_.size() - 1] / matrix_[matrix_.size() - 1][matrix_.size() - 1]);
		for (int i = matrix_.size() - 2; i >= 0; i--) {
			answers[i] = R((f_[i] - matrix_[i][i + 1] * answers[i + 1]) / matrix_[i][i]);
		}
		LogVector(answers);

		return answers;
	}

private:
	std::vector<std::vector<float>> matrix_;
	std::vector<float> f_, y_;
};

int main() {
	remove("logs.txt");
	Matrix m;
	std::cout << std::endl;
	std::vector<float> answers = m.Solve();
	std::cout << std::endl << "Vector of answers:" << std::endl;
	PrintVector(answers);
	std::cout << "Error is: " << CountInfelicity(m.GetY(), answers) << std::endl;
	std::cout << std::endl << std::endl;
	return 0;
}