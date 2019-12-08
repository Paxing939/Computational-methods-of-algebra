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
	Matrix(double E, int k_max) : E_(E), k_max_(k_max) {
		srand(time(NULL));
		// случайное количество элементов матрицы
		int n = rand() % 3 + 10;

		matrix_.resize(n);
		for (auto& line : matrix_) {
			line.resize(n);
		}

		// заполняем матрицу А
		// заполнение недиагональных элементов
		int p = 0, sum_1_col = 0;
		for (int i = 0; i < matrix_.size(); i++) {
			for (int j = 0; j < matrix_.size(); j++) {
				if (i != j) {
					matrix_[i][j] = rand() % 5 - 4;
				}
			}
		}

		// заполнение диагональных элементов
		for (int i = 1; i < matrix_.size(); i++) {
			matrix_[i][i] = 0;
			for (int j = 0; j < matrix_.size(); j++) {
				if (j != i) {
					matrix_[i][i] += -matrix_[i][j];
				}
			}
			sum_1_col += matrix_[0][i];
		}
		matrix_[0][0] = sum_1_col + 1;

		int m = 17;
		for (int i = 0; i < n; i++) {
			x_.push_back(m + i);
		}
		// с помощью умножения матриц находим вектора f
		f_ = MultMat(x_, matrix_);
		// пишем исходные матрицы в файл
		LogMatrix(matrix_);
	}

	bool EndAlgorithm(std::vector<double> x_new) {
		double max = 0.0f;
		for (int i = 0; i < x_new.size(); i++) {
			if (fabs(x_curr_[i] - x_new[i]) > max) {
				max = fabs(x_curr_[i] - x_new[i]);
			}
		}
		return max < E_ || k_ > k_max_;
	}

	std::vector<double> SolveJacobi() {
		// вывод на экран и в файл
		LogMatrix(matrix_); LogVector(f_);
		std::cout << "Vector x we need to get: " << std::endl;
		PrintVector(x_);
		std::cout << "Input matrix: " << std::endl;
		PrintMatrix(matrix_);
		std::cout << "Vector f: " << std::endl;
		PrintVector(f_);

		// метод Якоби
		x_curr_ = f_;
		std::vector<double> x_new(matrix_.size());
		for (k_ = 0; k_ < k_max_; k_++) {
			for (int i = 0; i < matrix_.size(); i++) {
				double sum = 0.0f; // считаем суммы для вычисления следующего x
				for (int j = 0; j < matrix_.size(); j++) {
					if (j != i) {
						sum += matrix_[i][j] * x_curr_[j];
					}
				}
				x_new[i] = 1 / matrix_[i][i] * (f_[i] - sum);
			}

			// если норма разности векторов меньше эпсилон - выходим из цикла
			if (EndAlgorithm(x_new)) {
				break;
			}
			x_curr_ = x_new;
		}
		// выводим сообщение если метод вышел по ограничению количества итераций
		if (k_ == k_max_) {
			std::cout << "Exceeding the maximum number of iterations!" << std::endl;
		}
		else {
			std::cout << "Number of iterations needed to success: " << k_ << std::endl;
		}

		return x_new;
	}

	std::vector<double> SolveGaussSeidel(double omega) {
		// вывод на экран и в файл
		LogMatrix(matrix_); LogVector(f_);
		std::cout << "Vector x we need to get: " << std::endl;
		PrintVector(x_);
		std::cout << "Input matrix: " << std::endl;
		PrintMatrix(matrix_);
		std::cout << "Vector f: " << std::endl;
		PrintVector(f_);

		// метод Гаусса-Зейделя
		x_curr_ = f_;
		std::vector<double> x_new = x_curr_;
		for (k_ = 0; k_ < k_max_; k_++) {
			for (int i = 0; i < matrix_.size(); i++) {
				double sum = 0.0f; // считаем суммы для вычисления следующего x
				for (int j = 0; j < i; j++) {
					sum += matrix_[i][j] * x_new[j];
				}
				for (int j = i + 1; j < matrix_.size(); j++) {
					sum += matrix_[i][j] * x_curr_[j];
				}
				x_new[i] = (1 - omega) * x_curr_[i] + omega / matrix_[i][i] * (f_[i] - sum);
			}
			// если норма разности векторов меньше эпсилон - выходим из цикла
			if (EndAlgorithm(x_new)) {
				break;
			}
			x_curr_ = x_new;
		}
		// выводим сообщение если метод вышел по ограничению количества итераций
		if (k_ == k_max_) {
			std::cout << "Exceeding the maximum number of iterations!" << std::endl;
		} else {
			std::cout << "Number of iterations needed to success: " << k_ << std::endl;
		}
		return x_new;
	}

private:
	std::vector<std::vector<double>> matrix_;
	std::vector<double> f_, x_, x_curr_;
	int k_max_, k_;
	double E_;
};

int main() {
	remove("logs.txt");
	Matrix m(0.0001, 1000);
	std::cout << std::endl;
	std::vector<double> answers = m.SolveJacobi();
	std::cout << "Vector of answers:" << std::endl;
	PrintVector(answers);

	std::cout << std::endl << std::endl;
	double omega = 0.5f;
	answers = m.SolveGaussSeidel(omega);
	std::cout << std::endl << "Vector of answers:" << std::endl;
	PrintVector(answers);

	std::cout << std::endl << std::endl;
	omega = 1.0f;
	answers = m.SolveGaussSeidel(omega);
	std::cout << "Vector of answers:" << std::endl;
	PrintVector(answers);

	std::cout << std::endl << std::endl;
	omega = 1.5f;
	answers = m.SolveGaussSeidel(omega);
	std::cout << "Vector of answers:" << std::endl;
	PrintVector(answers);

	return 0;
}