#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>

// 被積分函數 f(x) = x^2 * log(x)
double f(double x) {
    return x * x * std::log(x);
}

// 高斯積分 (使用給定的節點與權重)
double gaussian_quadrature(double a, double b, int n) {
    std::vector<double> x, w;

    // 根據 n 選擇節點與權重（對應到 [-1, 1] 的 Legendre 多項式）
    if (n == 3) {
        x = { -std::sqrt(3.0 / 5), 0.0, std::sqrt(3.0 / 5) };
        w = { 5.0 / 9, 8.0 / 9, 5.0 / 9 };
    } else if (n == 4) {
        x = {
            -std::sqrt((3 + 2 * std::sqrt(6.0 / 5)) / 7),
            -std::sqrt((3 - 2 * std::sqrt(6.0 / 5)) / 7),
             std::sqrt((3 - 2 * std::sqrt(6.0 / 5)) / 7),
             std::sqrt((3 + 2 * std::sqrt(6.0 / 5)) / 7)
        };
        w = {
            (18 - std::sqrt(30)) / 36,
            (18 + std::sqrt(30)) / 36,
            (18 + std::sqrt(30)) / 36,
            (18 - std::sqrt(30)) / 36
        };
    } else {
        throw std::invalid_argument("Only n = 3 or 4 supported in this implementation.");
    }

    // 範圍變換 [a, b]
    double integral = 0.0;
    for (int i = 0; i < n; ++i) {
        double t = 0.5 * (x[i] + 1) * (b - a) + a;
        integral += w[i] * f(t);
    }

    return 0.5 * (b - a) * integral;
}

// 真值（手動計算或使用其他方法求得）
double exact_integral() {
    // f(x) = x^2 * log(x) 的不定積分為：
    // (x^3 / 3) * log(x) - x^3 / 9
    auto F = [](double x) {
        return (x * x * x / 3.0) * std::log(x) - x * x * x / 9.0;
    };
    return F(1.5) - F(1.0);
}

int main() {
    double a = 1.0;
    double b = 1.5;

    double exact = exact_integral();
    double gauss3 = gaussian_quadrature(a, b, 3);
    double gauss4 = gaussian_quadrature(a, b, 4);

    std::cout << std::fixed << std::setprecision(10);
    std::cout << "Exact value (analytical):          " << exact << std::endl;
    std::cout << "Gaussian Quadrature (n=3):         " << gauss3 << std::endl;
    std::cout << "Gaussian Quadrature (n=4):         " << gauss4 << std::endl;
    std::cout << "Error (n=3):                       " << std::abs(gauss3 - exact) << std::endl;
    std::cout << "Error (n=4):                       " << std::abs(gauss4 - exact) << std::endl;

    return 0;
}
