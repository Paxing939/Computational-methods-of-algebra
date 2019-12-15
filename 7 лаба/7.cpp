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

std::vector<float> DecreaseCoeffs(std::vector<float> coeffs) {
	for (auto& coeff : coeffs) {
		coeff /= coeffs[coeffs.size() - 1];
	}
	return coeffs;
}

std::vector<float> Solve2Equation(std::vector<float> coeffs) {
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

bool OneSign(float first, float second) {
	if (first < 0 && second < 0) {
		return true;
	}
	if (first > 0 && second < 0) {
		return false;
	}
	if (first < 0 && second > 0) {
		return false;
	}
	if (first > 0 && second > 0) {
		return true;
	}
}

//std::vector<float> SeparateRoots(std::vector<float> coeffs, float extr) {
//	float point_2 = 0, margin = 10;
//	for (int i = 0; i < 100; i++) {
//		point_2 = CountPolynomFromX(coeffs, extr - margin);
//		if (OneSign(extr, point_2)) {
//			margin++;
//			continue;
//		}
//		if (!OneSign(extr - 1 , point_2)) {
//			extr--;
//		}
//		else {
//			break;
//		}
//	}
//}

float XFromPolynom(std::vector<float> coeffs, float x) {
	return -(coeffs[3] * x * x * x + coeffs[2] * x * x + coeffs[0]) / coeffs[1];
}

float XFromPolynom2(std::vector<float> coeffs, float x) {
	//return -(coeffs[3] * x * x * x + coeffs[2] * x * x + coeffs[0]) / coeffs[1];
	return -pow((coeffs[2] * x * x + coeffs[1] * x + coeffs[0]) / coeffs[3], 1 / 3);
	return (coeffs[2] + coeffs[1] / x + coeffs[0] / (x * x)) / coeffs[3];
}

float XFromPolynom22(std::vector<float> coeffs, float x) {
	//return -(coeffs[3] * x * x * x + coeffs[2] * x * x + coeffs[0]) / coeffs[1];
	//return -pow((coeffs[2] * x * x + coeffs[1] * x + coeffs[0]) / coeffs[3], 1 / 3);
	return (coeffs[2] + coeffs[1] / x + coeffs[1] / (x * x)) / coeffs[3];
}

float XFromPolynom3(std::vector<float> coeffs, float x) {
	return (coeffs[2] + coeffs[1] / x + coeffs[1] / (x * x)) / coeffs[3];
}

class Polinomial {
public:

	~Polinomial() {
		coeffs_.clear();
	}

	Polinomial() {
		coeffs_ = { -5, -20, -4, 4, 1 };
		accurancy_ = 0.0000001;
	}

	void StartIterations() {
		std::vector<float> derivative_3 = DecreaseCoeffs(FindDerivative(coeffs_));
		std::vector<float> derivative_2 = FindDerivative(derivative_3);
		std::vector<float> derivative_2_answers = Solve2Equation(derivative_2), derivative_3_answers, answers;
		if (derivative_2_answers[0] > derivative_2_answers[1]) {
			std::swap(derivative_2_answers[0], derivative_2_answers[1]);
		}

		//#define DEBAG 1
		//		derivative_3 = { 5, -8, -8, 5 };
		//		derivative_2_answers = { -0.37098, 1.4376 };

		float margin = 5;
		float x_2 = derivative_2_answers[0], x_1 = derivative_2_answers[0] - margin;
		float y_1 = CountPolynomFromX(derivative_3, x_1), y_2 = CountPolynomFromX(derivative_3, x_2);
		for (int i = 0; i < 100; i++) {
			if (OneSign(y_1, y_2)) {
				x_2--;
				y_2 = CountPolynomFromX(derivative_3, x_2);
			}
			else {
				break;
			}
		}
		float x = (x_1 + x_2) / 2, x_next = 0;
		for (int i = 0; i < 1500; i++) {
			x_next = XFromPolynom2(derivative_3, x);
			if (abs(x_next - x) < accurancy_) {
				break;
			}
			x = x_next;
		}
		std::cout << x << std::endl;
		derivative_3_answers.push_back(x);
		std::cout << x << std::endl;

		x_1 = derivative_2_answers[0], x_2 = derivative_2_answers[1];
		y_1 = CountPolynomFromX(derivative_3, x_1), y_2 = CountPolynomFromX(derivative_3, x_2);
		if (OneSign(y_1, y_2)) {
			std::cout << "Complex root in f'(x)!" << std::endl;
			exit(0);
		}
		//for (int i = 0; i < 100; i++) {
		//	if (!OneSign(CountPolynomFromX(derivative_3, x_1 + 1), y_2)) {
		//		x_1++;
		//		y_1 = CountPolynomFromX(derivative_3, x_1);
		//	} else {
		//		if (!OneSign(CountPolynomFromX(derivative_3, x_2 - 1), y_1)) {
		//			x_2--;
		//			y_1 = CountPolynomFromX(derivative_3, x_2);
		//		} else {
		//			break;
		//		}
		//	}
		//}
		x = (x_2 - x_1) / 2, x_next = 0;
		for (int i = 0; i < 500; i++) {
			x_next = XFromPolynom2(derivative_3, x);
			if (abs(x_next - x) < accurancy_) {
				break;
			}
			x = x_next;
		}
		derivative_3_answers.push_back(x);
		std::cout << x << std::endl;


		x_1 = derivative_2_answers[1], x_2 = derivative_2_answers[1] + 10;
		y_1 = CountPolynomFromX(derivative_3, x_1), y_2 = CountPolynomFromX(derivative_3, x_2);
		for (int i = 0; i < 100; i++) {
			if (OneSign(y_1, y_2)) {
				x_2++;
				y_2 = CountPolynomFromX(derivative_3, x_2);
			}
		}
		x = x_1, x_next = 0;
		for (int i = 0; i < 500; i++) {
			x_next = XFromPolynom2(derivative_3, x);
			if (abs(x_next - x) < accurancy_) {
				break;
			}
			x = x_next;
		}
		derivative_3_answers.push_back(x);
		std::cout << x << std::endl << std::endl;




		x_2 = derivative_3_answers[0], x_1 = derivative_3_answers[0] - margin;
		y_1 = CountPolynomFromX(coeffs_, x_1), y_2 = CountPolynomFromX(coeffs_, x_2);
		for (int i = 0; i < 100; i++) {
			if (OneSign(y_1, y_2)) {
				x_1--;
			}
			y_1 = CountPolynomFromX(coeffs_, x_1);
		}
		x = (x_1 - x_2) / 2, x_next = 0;
		for (int i = 0; i < 1500; i++) {
			x_next = XFromPolynom22(coeffs_, x);
			if (abs(x_next - x) < accurancy_) {
				break;
			}
			x = x_next;
		}
		answers.push_back(x);
		std::cout << x << std::endl;

		x_2 = derivative_3_answers[0], x_1 = derivative_3_answers[1];
		y_1 = CountPolynomFromX(coeffs_, x_1), y_2 = CountPolynomFromX(coeffs_, x_2);
		x = (x_1 - x_2) / 2, x_next = 0;
		for (int i = 0; i < 1500; i++) {
			x_next = XFromPolynom22(coeffs_, x);
			if (abs(x_next - x) < accurancy_) {
				break;
			}
			x = x_next;
		}
		answers.push_back(x);
		std::cout << x << std::endl;

		x_2 = derivative_3_answers[1], x_1 = derivative_3_answers[2];
		y_1 = CountPolynomFromX(coeffs_, x_1), y_2 = CountPolynomFromX(coeffs_, x_2);
		x = (x_1 - x_2) / 2, x_next = 0;
		for (int i = 0; i < 1500; i++) {
			x_next = XFromPolynom22(coeffs_, x);
			if (abs(x_next - x) < accurancy_) {
				break;
			}
			x = x_next;
		}
		std::cout << x << std::endl;
		answers.push_back(x);
		std::cout << x << std::endl;

		x_2 = derivative_3_answers[2], x_1 = derivative_3_answers[2] + margin;
		y_1 = CountPolynomFromX(coeffs_, x_1), y_2 = CountPolynomFromX(coeffs_, x_2);
		for (int i = 0; i < 100; i++) {
			if (OneSign(y_1, y_2)) {
				x_1++;
			}
			y_1 = CountPolynomFromX(coeffs_, x_1);
		}
		x = (x_1 - x_2) / 2, x_next = 0;
		for (int i = 0; i < 1500; i++) {
			x_next = XFromPolynom22(coeffs_, x);
			if (abs(x_next - x) < accurancy_) {
				break;
			}
			x = x_next;
		}
		std::cout << x << std::endl;
		answers.push_back(x);
		std::cout << x << std::endl << std::endl << std::endl << std::endl;
	}

	void Newton() {
		//std::vector<float> derivative = DecreaseCoeffs(FindDerivative(coeffs_));
		//float x = 1, x_next = 0;
		//for (int i = 0; i < 100; i++) {
		//	x_next = x - CountPolynomFromX(coeffs_, x) / CountPolynomFromX(derivative, x);
		//	if (abs(x_next - x) < accurancy_) {
		//		break;
		//	}
		//}


		std::vector<float> derivative_3 = DecreaseCoeffs(FindDerivative(coeffs_));
		std::vector<float> derivative_2 = FindDerivative(derivative_3);
		std::vector<float> derivative_2_answers = Solve2Equation(derivative_2), derivative_3_answers, answers;
		if (derivative_2_answers[0] > derivative_2_answers[1]) {
			std::swap(derivative_2_answers[0], derivative_2_answers[1]);
		}

		//#define DEBAG 1
		//		derivative_3 = { 5, -8, -8, 5 };
		//		derivative_2_answers = { -0.37098, 1.4376 };

		float margin = 5;
		float x_2 = derivative_2_answers[0], x_1 = derivative_2_answers[0] - margin;
		float y_1 = CountPolynomFromX(derivative_3, x_1), y_2 = CountPolynomFromX(derivative_3, x_2);
		for (int i = 0; i < 100; i++) {
			if (OneSign(y_1, y_2)) {
				x_2--;
				y_2 = CountPolynomFromX(derivative_3, x_2);
			}
			else {
				break;
			}
		}
		float x = x_1, x_next = 0;
		for (int i = 0; i < 1500; i++) {
			x_next = x - CountPolynomFromX(coeffs_, x) / CountPolynomFromX(derivative_3, x);
			if (abs(x_next - x) < accurancy_) {
				break;
			}
			x = x_next;
		}
		std::cout << x << std::endl;
		derivative_3_answers.push_back(x);

		x_1 = derivative_2_answers[0], x_2 = derivative_2_answers[1];
		y_1 = CountPolynomFromX(derivative_3, x_1), y_2 = CountPolynomFromX(derivative_3, x_2);
		if (OneSign(y_1, y_2)) {
			std::cout << "Complex root in f'(x)!" << std::endl;
			exit(0);
		}
		//for (int i = 0; i < 100; i++) {
		//	if (!OneSign(CountPolynomFromX(derivative_3, x_1 + 1), y_2)) {
		//		x_1++;
		//		y_1 = CountPolynomFromX(derivative_3, x_1);
		//	} else {
		//		if (!OneSign(CountPolynomFromX(derivative_3, x_2 - 1), y_1)) {
		//			x_2--;
		//			y_1 = CountPolynomFromX(derivative_3, x_2);
		//		} else {
		//			break;
		//		}
		//	}
		//}
		x = x_1, x_next = 0;
		for (int i = 0; i < 500; i++) {
			x_next = x - CountPolynomFromX(coeffs_, x) / CountPolynomFromX(derivative_3, x);
			if (abs(x_next - x) < accurancy_) {
				break;
			}
			x = x_next;
		}
		derivative_3_answers.push_back(x);
		std::cout << x << std::endl;


		x_1 = derivative_2_answers[1], x_2 = derivative_2_answers[1] + 10;
		y_1 = CountPolynomFromX(derivative_3, x_1), y_2 = CountPolynomFromX(derivative_3, x_2);
		for (int i = 0; i < 100; i++) {
			if (OneSign(y_1, y_2)) {
				x_2++;
				y_2 = CountPolynomFromX(derivative_3, x_2);
			}
		}
		x = x_1, x_next = 0;
		for (int i = 0; i < 500; i++) {
			x_next = x - CountPolynomFromX(coeffs_, x) / CountPolynomFromX(derivative_3, x);
			if (abs(x_next - x) < accurancy_) {
				break;
			}
			x = x_next;
		}
		derivative_3_answers.push_back(x);
		std::cout << x << std::endl << std::endl;




		x_2 = derivative_3_answers[0], x_1 = derivative_3_answers[0] - margin;
		y_1 = CountPolynomFromX(coeffs_, x_1), y_2 = CountPolynomFromX(coeffs_, x_2);
		for (int i = 0; i < 100; i++) {
			if (OneSign(y_1, y_2)) {
				x_1--;
			}
			y_1 = CountPolynomFromX(coeffs_, x_1);
		}
		x = (x_1 - x_2) / 2, x_next = 0;
		for (int i = 0; i < 1500; i++) {
			x_next = x - CountPolynomFromX(coeffs_, x) / CountPolynomFromX(derivative_3, x);
			if (abs(x_next - x) < accurancy_) {
				break;
			}
			x = x_next;
		}
		answers.push_back(x);
		std::cout << x << std::endl;

		x_2 = derivative_3_answers[0], x_1 = derivative_3_answers[1];
		y_1 = CountPolynomFromX(coeffs_, x_1), y_2 = CountPolynomFromX(coeffs_, x_2);
		x = x_1, x_next = 0;
		for (int i = 0; i < 1500; i++) {
			x_next = x - CountPolynomFromX(coeffs_, x) / CountPolynomFromX(derivative_3, x);
			if (abs(x_next - x) < accurancy_) {
				break;
			}
			x = x_next;
		}
		std::cout << x << std::endl;
		answers.push_back(x);

		x_2 = derivative_3_answers[1], x_1 = derivative_3_answers[2];
		y_1 = CountPolynomFromX(coeffs_, x_1), y_2 = CountPolynomFromX(coeffs_, x_2);
		x = x_1, x_next = 0;
		for (int i = 0; i < 1500; i++) {
			x_next = x - CountPolynomFromX(coeffs_, x) / CountPolynomFromX(derivative_3, x);
			if (abs(x_next - x) < accurancy_) {
				break;
			}
			x = x_next;
		}
		std::cout << x << std::endl;
		answers.push_back(x);

		x_2 = derivative_3_answers[2], x_1 = derivative_3_answers[2] + margin;
		y_1 = CountPolynomFromX(coeffs_, x_1), y_2 = CountPolynomFromX(coeffs_, x_2);
		for (int i = 0; i < 100; i++) {
			if (OneSign(y_1, y_2)) {
				x_1++;
			}
			y_1 = CountPolynomFromX(coeffs_, x_1);
		}
		x = x_1, x_next = 0;
		for (int i = 0; i < 1500; i++) {
			x_next = x - CountPolynomFromX(coeffs_, x) / CountPolynomFromX(derivative_3, x);
			if (abs(x_next - x) < accurancy_) {
				break;
			}
			x = x_next;
		}
		std::cout << x << std::endl;
		answers.push_back(x);
	}

private:
	float accurancy_;
	std::vector<float> coeffs_;
};

int main() {
	remove("logs.txt");
	Polinomial m;
	std::cout << std::endl;
	
	m.Newton();
	m.StartIterations();
	return 0;
}