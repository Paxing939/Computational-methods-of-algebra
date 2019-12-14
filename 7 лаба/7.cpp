#include <cmath>
#include <iomanip>
#include "Functions.h"

std::vector<float> FindDerivative(std::vector<float> coeffs) {
	std::vector<float> derivative;
	for (int i = 1; i < coeffs.size(); i++) {
		derivative.push_back(coeffs[i] * i);
	}
	return derivative;
}

std::pair<float, float> Solve2Equation(std::vector<float> coeffs) {
	float D = coeffs[1] * coeffs[1] - 4 * coeffs[2] * coeffs[0];
	if (D < 0) {
		std::cout << "Discriminat < 0!" << std::endl;
		exit(0);
	}
	return { (-coeffs[1] + sqrt(D)) / (2 * coeffs[2]), (-coeffs[1] - sqrt(D)) / (2 * coeffs[2]) };
}


float CountPolynomFromX(std::vector<float> coeffs, float x) {
	float res = 0;
	for (int i = 0; i < coeffs.size(); i++) {
		res += coeffs[i] * pow(x, i);
	}
	return res;
}

class Polinomial {
public:

	~Polinomial() {
		coeffs_.clear();
	}

	// создание матрицы А
	Polinomial() {
		coeffs_ = { -8, 9, 3, -5, 12 };
	}

	void StartIterations() {
		std::vector<float> derivative_3 = FindDerivative(coeffs_);
		std::vector<float> derivative_2 = FindDerivative(coeffs_);
		std::pair<float, float> derivative_2_answers = Solve2Equation(derivative_2);
		for (int i = 0; i < 3; i++) {

		}

	}

private:
	std::vector<float> coeffs_;
};

int main() {
	remove("logs.txt");
	Polinomial m;
	std::cout << std::endl;
	m.StartIterations();

	return 0;
}