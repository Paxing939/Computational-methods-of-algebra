#include <vector>
#include <iostream>
#include <fstream>

double R(double x) {
    int k = 1000;
    return round(x * k) / (double)k;
}

double CountInfelicity(std::vector<double> x_de_ure, std::vector<double> x_de_facto) {
    double sum_x2 = 0, sum_x = 0;
    for (int i = 0; i < x_de_facto.size(); i++) {
        sum_x2 += x_de_facto[i] * x_de_facto[i];
        sum_x += (x_de_facto[i] - x_de_ure[i]) * (x_de_facto[i] - x_de_ure[i]);
    }
    std::ofstream fout("errors.txt");
    fout << sum_x / sum_x2 << std::endl;
    fout.close();
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


// функция умножения матрицы на вектор
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

void PrintMatrix(const std::vector<std::vector<double>>& matrix) {
    for (const auto& vec : matrix) {
        for (const auto& x_ : vec) {
            std::cout << std::setw(7) << std::setprecision(2) << std::left << x_ << " ";
        }
        std::cout << std::endl;
    }
}

void PrintVector(const std::vector<double>& vector) {
    for (const auto& x_ : vector) {
        std::cout << std::setw(7) << std::setprecision(2) << std::left << x_ << " ";
    }
    std::cout << std::endl;
}
