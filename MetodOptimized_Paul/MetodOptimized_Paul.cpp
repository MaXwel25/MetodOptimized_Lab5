#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>

using namespace std;

// главная функция f(x) = x1^2 + 6*x2^2 + x1*x2 + x1
long double main_function(double x1, double x2) {
    return x1 * x1 + 6 * x2 * x2 + x1 * x2 + x1;
}

vector<double> gradient(double x1, double x2) {
    return { 2 * x1 + x2 + 1, 12 * x2 + x1 };
}

// норма вектора
long double norma(const vector<double>& v) {
    return sqrt(v[0] * v[0] + v[1] * v[1]);
}

// умножение матрицы 2x2 на вектор
vector<double> mat_vec_mult(const vector<vector<double>>& A, const vector<double>& v) {
    return {
        A[0][0] * v[0] + A[0][1] * v[1],
        A[1][0] * v[0] + A[1][1] * v[1]
    };
}

// внешнее произведение векторов
vector<vector<double>> outer_product(const vector<double>& u, const vector<double>& v) {
    return {
        {u[0] * v[0], u[0] * v[1]},
        {u[1] * v[0], u[1] * v[1]}
    };
}

// сложение матриц
vector<vector<double>> matrix_add(const vector<vector<double>>& A, const vector<vector<double>>& B) {
    return {
        {A[0][0] + B[0][0], A[0][1] + B[0][1]},
        {A[1][0] + B[1][0], A[1][1] + B[1][1]}
    };
}

// вычитание матриц
vector<vector<double>> matrix_sub(const vector<vector<double>>& A, const vector<vector<double>>& B) {
    return {
        {A[0][0] - B[0][0], A[0][1] - B[0][1]},
        {A[1][0] - B[1][0], A[1][1] - B[1][1]}
    };
}

// умножение матрицы на скаляр
vector<vector<double>> matrix_scale(const vector<vector<double>>& A, double s) {
    return {
        {A[0][0] * s, A[0][1] * s},
        {A[1][0] * s, A[1][1] * s}
    };
}

// скалярное произведение
double skalar(const vector<double>& u, const vector<double>& v) {
    return u[0] * v[0] + u[1] * v[1];
}

// метод золотого сечения для одномерной минимизации (для 11 шага)
double golden_section(double x1, double x2, const vector<double>& direction,
    double a, double b, double eps = 1e-8) {
    const double phi = (1.0 + sqrt(5.0)) / 2.0;
    double x1_val, x2_val;

    while (fabs(b - a) > eps) {
        x1_val = b - (b - a) / phi;
        x2_val = a + (b - a) / phi;
        double f1 = main_function(x1 + x1_val * direction[0], x2 + x1_val * direction[1]);
        double f2 = main_function(x1 + x2_val * direction[0], x2 + x2_val * direction[1]);
        if (f1 > f2) {
            a = x1_val;
        }
        else {
            b = x2_val;
        }
    }

    return (a + b) / 2.0;
}

// метод Дэвидона-Флетчера-Пауэлла
vector<double> dfp_method(double x1_0, double x2_0, double eps1 = 0.1, double eps2 = 0.15, int M = 10) {
    // 1/2) инициализация
    vector<double> x = { x1_0, x2_0 };
    vector<vector<double>> A = { {1.0, 0.0}, {0.0, 1.0} }; // A^0 = E
    int k = 0;
    vector<double> x_prev_iter = x;
    vector<double> grad_prev;
    vector<double> x_prev = x;
    double f_prev = main_function(x[0], x[1]);
    cout << "Метод Дэвидона-Флетчера-Пауэлла\n";
    cout << "Функция: f(x1, x2) = x1^2 + 6*x2^2 + x1*x2 + x1" << endl;
    cout << "Начальная точка: (" << x[0] << ", " << x[1] << ")" << endl;
    cout << "f(x0) = " << f_prev << endl;
    cout << "Параметры: e1 = " << eps1 << ", e2 = " << eps2 << ", M = " << M << "\n\n";

    while (true) {
        // 3) вычислить градиент
        vector<double> grad = gradient(x[0], x[1]);

        // 4) проверяем критерий окончания
        double grad_norm = norma(grad);
        cout << "Итерация " << k << ":\n";
        cout << "x = (" << setprecision(6) << x[0] << ", " << x[1] << ")\n";
        cout << "f(x) = " << setprecision(6) << main_function(x[0], x[1]) << endl;
        cout << "*f(x) = (" << setprecision(6) << grad[0] << ", " << grad[1] << ")\n";
        cout << "||*f(x)|| = " << setprecision(6) << grad_norm << endl;

        if (grad_norm < eps1) {
            cout << "Критерий окончания ||*f(x)|| < e1 выполнен" << endl;
            break;
        }
        // 5) проверяем максимальное число итераций
        if (k >= M) {
            cout << "Достигнуто максимальное число итераций" << endl;
            break;
        }

        vector<double> main_vec;

        if (k == 0) {
            // при k = 0 перейти к шагу 10
            // 10) определяем направление
            main_vec = mat_vec_mult(A, grad);
            main_vec[0] = -main_vec[0];
            main_vec[1] = -main_vec[1];
        }
        else {
            // 6) вычислить *g^k
            vector<double> delta_g = { grad[0] - grad_prev[0], grad[1] - grad_prev[1] };

            // 7) вычислить *x^k
            vector<double> delta_x = { x[0] - x_prev_iter[0], x[1] - x_prev_iter[1] };

            // 8) вычислить A^k_c
            double delta_x_dot_delta_g = skalar(delta_x, delta_g);
            vector<double> A_delta_g = mat_vec_mult(A, delta_g);
            double delta_g_dot_A_delta_g = skalar(delta_g, A_delta_g);

            if (fabs(delta_x_dot_delta_g) > 1e-15 && fabs(delta_g_dot_A_delta_g) > 1e-15) {
                vector<vector<double>> term1 = outer_product(delta_x, delta_x);
                term1 = matrix_scale(term1, 1.0 / delta_x_dot_delta_g);

                vector<vector<double>> term2 = outer_product(A_delta_g, A_delta_g);
                term2 = matrix_scale(term2, 1.0 / delta_g_dot_A_delta_g);

                vector<vector<double>> A_c = matrix_sub(term1, term2);

                // 9) обновить матрицу A
                A = matrix_add(A, A_c);
            }

            // 10) определить направление
            main_vec = mat_vec_mult(A, grad);
            main_vec[0] = -main_vec[0];
            main_vec[1] = -main_vec[1];
        }

        // 11) одномерная минимизация
        double t = golden_section(x[0], x[1], main_vec, 0.0, 1.0);

        // сохраняем значения для следующей итерации
        x_prev_iter = x;
        grad_prev = grad;

        // 12) вычислить новую точку
        x_prev = x;
        f_prev = main_function(x[0], x[1]);
        x[0] = x[0] + t * main_vec[0];
        x[1] = x[1] + t * main_vec[1];

        double f_current = main_function(x[0], x[1]);

        cout << "t = " << t << endl;
        cout << "Направление: (" << main_vec[0] << ", " << main_vec[1] << ")\n";
        cout << "Матрица A:\n";
        cout << "[" << A[0][0] << "  " << A[0][1] << "]\n";
        cout << "[" << A[1][0] << "  " << A[1][1] << "]\n";

        // 13) проверить дополнительные критерии окончания
        if (k >= 1) {
            double dx_norm = norma({ x[0] - x_prev[0], x[1] - x_prev[1] });
            double df_norm = fabs(f_current - f_prev);

            cout << "||*x|| = " << dx_norm;
            cout << ", |*f| = " << df_norm << endl;

            if (dx_norm < eps2 && df_norm < eps2) {
                cout << "Критерии ||x^(k+1) - x^k|| < e2 и |f(x^(k+1)) - f(x^k)| < e2 выполнены\n";
                break;
            }
        }
        cout << "----------------------------------------\n";
        k++;
    }

    return x;
}

int main() {
    setlocale(LC_ALL, "ru");

    // параметры метода
    long double x1_0 = 1.5, x2_0 = 1.1;
    long double eps1 = 0.1, eps2 = 0.15;
    int M = 10;

    vector<double> solution = dfp_method(x1_0, x2_0, eps1, eps2, M);

    // вывод результатов
    cout << "\nИтоговые результаты:\n";
    cout << "Точка минимума: (" << solution[0]
        << ", " << solution[1] << ")\n";
    cout << "Значение функции: " << main_function(solution[0], solution[1]) << endl;

    vector<double> final_grad = gradient(solution[0], solution[1]);
    cout << "Градиент: (" << final_grad[0]
        << ", " << final_grad[1] << ")" << endl;
    cout << "Норма градиента: " << norma(final_grad) << endl;

    // аналитическое решение для проверки
    long double x1_analytical = -12.0 / 23.0;
    long double x2_analytical = 1.0 / 23.0;
    cout << endl << "Аналитическое решение:" << endl;
    cout << "x1* = " << x1_analytical << endl;
    cout << "x2* = " << x2_analytical << endl;
    cout << "f(x*) = " << main_function(x1_analytical, x2_analytical) << endl;

    return 0;
}