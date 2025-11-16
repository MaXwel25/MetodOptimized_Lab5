#pragma once

#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

// целевая функция: f(x) = x1^2 + 6x2^2 + x1x2 + x1
double main_funct(const vector<double>& x) {
    return x[0] * x[0] + 6 * x[1] * x[1] + x[0] * x[1] + x[0];
}

// градиент функции
vector<double> gradient(const vector<double>& x) {
    vector<double> grad(2);
    grad[0] = 2 * x[0] + x[1] + 1;
    grad[1] = x[0] + 12 * x[1];
    return grad;
}

// скалярное произведение векторов
double skalar_proiz(const vector<double>& a, const vector<double>& b) {
    double result = 0.0;
    for (size_t i = 0; i < a.size(); i++) {
        result += a[i] * b[i];
    }
    return result;
}

// норма вектора
double vector_norm(const vector<double>& v) {
    return sqrt(skalar_proiz(v, v));
}

// умножение вектора на скаляр
vector<double> scalar_multiply(const vector<double>& v, double scalar) {
    vector<double> result = v;
    for (double& val : result) {
        val *= scalar;
    }
    return result;
}

// сложение векторов
vector<double> vector_add(const vector<double>& a, const vector<double>& b) {
    vector<double> result = a;
    for (size_t i = 0; i < a.size(); i++) {
        result[i] += b[i];
    }
    return result;
}

// вывод вектора
void print_vector_pollak(const vector<double>& v, const string& name = "") {
    if (!name.empty()) {
        cout << name << " = ";
    }
    cout << "(";
    for (size_t i = 0; i < v.size(); i++) {
        cout << v[i] << (i < v.size() - 1 ? ", " : "");
    }
    cout << ")";
}

// одномерный поиск методом золотого сечения
double golden_section(const vector<double>& x0,
    const vector<double>& direction,
    double a = -1.0, double b = 1.0) {
    const double golden_ratio = (sqrt(5.0) - 1.0) / 2.0;
    const double eps = 1e-8;
    int n = x0.size();

    double c = b - golden_ratio * (b - a);
    double d = a + golden_ratio * (b - a);

    auto point = [&](double lambda) {
        vector<double> x = x0;
        for (int i = 0; i < n; i++) {
            x[i] += lambda * direction[i];
        }
        return main_funct(x);
        };

    double fc = point(c);
    double fd = point(d);

    while (abs(b - a) > eps) {
        if (fc < fd) {
            b = d;
            d = c;
            fd = fc;
            c = b - golden_ratio * (b - a);
            fc = point(c);
        }
        else {
            a = c;
            c = d;
            fc = fd;
            d = a + golden_ratio * (b - a);
            fd = point(d);
        }
    }
    return (a + b) / 2.0;
}

// метод сопряженных градиентов (Поллака-Райбера)
vector<double> Pollaka(const vector<double>& x0,
    double tol = 1e-6,
    int max_iter = 1000) {
    vector<double> x = x0;
    int n = x.size();

    cout << "\n\nМетод Поллака-Райбера\n";
    cout << "Целевая функция: f(x) = x1^2 + 6x2^2 + x1x2 + x1\n";
    cout << "Начальная точка: ";
    print_vector_pollak(x0);
    cout << ", f = " << main_funct(x) << endl;
    cout << "Точность: " << tol << endl << endl;

    // начальный градиент и направление
    vector<double> g = gradient(x);
    vector<double> d = scalar_multiply(g, -1.0);
    double g_norm = vector_norm(g);

    cout << "Итерация 0:\n";
    cout << "  x0 = "; print_vector_pollak(x); cout << endl;
    cout << "  f(x0) = " << main_funct(x) << endl;
    cout << "  *f(x0) = "; print_vector_pollak(g); cout << endl;
    cout << "  |*f(x0)| = " << g_norm << endl;
    cout << "  d0 = "; print_vector_pollak(d); cout << endl << endl;

    for (int k = 0; k < max_iter; k++) {
        // одномерная минимизация вдоль направления d
        double lambda = golden_section(x, d, -2.0, 2.0);

        // сохраняем старые значения
        vector<double> x_old = x;
        vector<double> g_old = g;

        // обновляем точку
        x = vector_add(x, scalar_multiply(d, lambda));

        // вычисляем новый градиент
        g = gradient(x);
        g_norm = vector_norm(g);
        cout << "Итерация " << k + 1 << ":\n";
        cout << "  lambda = " << lambda << endl;
        cout << "  x = "; print_vector_pollak(x); cout << endl;
        cout << "  f(x) = " << main_funct(x) << endl;
        cout << "  *f(x) = "; print_vector_pollak(g); cout << endl;
        cout << "  |*f(x)| = " << g_norm << endl;

        // проверка критерия остановки
        if (g_norm < tol) {
            cout << "\nСходимость достигнута на итерации " << k + 1 << " !!!" << endl;
            break;
        }

        // вычисление коэффициента b (Флетчера-Ривса)
        double beta = skalar_proiz(g, g) / skalar_proiz(g_old, g_old);
        cout << "  b = " << beta << endl;

        // Обновление направления
        d = vector_add(scalar_multiply(g, -1.0), scalar_multiply(d, beta));
        cout << "  d = "; print_vector_pollak(d); cout << endl << endl;
        if (k == max_iter - 1) {
            cout << "\nДостигнуто максимальное число итераций\n";
        }
    }
    return x;
}