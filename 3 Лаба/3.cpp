#include <cmath>
#include <iomanip>
#include "Functions.h"

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

		f_ = MultMat(matrix_, y_);
	}

	std::vector<double> GetY() { return y_; }

	std::vector<double> Solve() {
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
			double dev = R(matrix_[i][i - 1] / matrix_[i - 1][i - 1]);
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
		std::vector<double> answers(matrix_.size());
		answers[answers.size() - 1] = R(f_[f_.size() - 1] / matrix_[matrix_.size() - 1][matrix_.size() - 1]);
		for (int i = matrix_.size() - 2; i >= 0; i--) {
			answers[i] = R((f_[i] - matrix_[i][i + 1] * answers[i + 1]) / matrix_[i][i]);
		}
		LogVector(answers);

		return answers;
	}

private:
	std::vector<std::vector<double>> matrix_;
	std::vector<double> f_, y_;
};

int main() {
	remove("logs.txt");
	Matrix m;
	std::cout << std::endl;
	std::vector<double> answers = m.Solve();
	std::cout << std::endl << "Vector of answers:" << std::endl;
	PrintVector(answers);
	std::cout << "Error is: " << CountInfelicity(m.GetY(), answers) << std::endl;
	std::cout << std::endl << std::endl;
	return 0;
}