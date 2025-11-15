#include <iostream>
#include <vector>
#include <cmath>
#include <functional>
using namespace std;

// целевая функция: f(x) = x1^2 + 6x2^2 + x1x2 + x1
double main_function(const vector<double>& x) {
    return x[0] * x[0] + 6 * x[1] * x[1] + x[0] * x[1] + x[0];
} 

// одномерный поиск методом золотого сечения
double golden_section(const function<double(const vector<double>&)>& func,
    const vector<double>& x0,
    const vector<double>& direction,
    double a = -1.0, double b = 1.0) {
    const double phi = (sqrt(5.0) - 1.0) / 2.0;
    const double eps = 1e-8;
    int n = x0.size();

    double c = b - phi * (b - a);
    double d = a + phi * (b - a);

    auto pointAlongDirection = [&](double lambda) {
        vector<double> x = x0;
        for (int i = 0; i < n; i++) {
            x[i] += lambda * direction[i];
        }
        return func(x);
        };

    double fc = pointAlongDirection(c);
    double fd = pointAlongDirection(d);

    while (abs(b - a) > eps) {
        if (fc < fd) {
            b = d;
            d = c;
            fd = fc;
            c = b - phi * (b - a);
            fc = pointAlongDirection(c);
        }
        else {
            a = c;
            c = d;
            fc = fd;
            d = a + phi * (b - a);
            fd = pointAlongDirection(d);
        }
    }
    return (a + b) / 2.0;
}

// нормализация вектора
vector<double> normal(const vector<double>& v) {
    int n = v.size();
    double norm = 0.0;
    for (double val : v) {
        norm += val * val;
    }
    norm = sqrt(norm);

    if (norm < 1e-12) return v;

    vector<double> result = v;
    for (double& val : result) {
        val /= norm;
    }
    return result;
}

// вычисление нормы вектора
double vector_norma(const vector<double>& v) {
    double norm = 0.0;
    for (double val : v) {
        norm += val * val;
    }
    return sqrt(norm);
}

// Печать вектора
void printVector(const vector<double>& v, const string& name = "") {
    if (!name.empty()) {
        cout << name << " = ";
    }
    cout << "(";
    for (size_t i = 0; i < v.size(); i++) {
        cout << v[i] << (i < v.size() - 1 ? ", " : "");
    }
    cout << ")";
}

// основная функция метода Пауэлла
vector<double> Paul(const function<double(const vector<double>&)>& func,
    const vector<double>& x0,
    double tol = 1e-6,
    int max_iter = 1000) {
    int n = x0.size();
    vector<double> x = x0;

    // инициализация направлений (координатные оси)
    vector<vector<double>> directions;
    for (int i = 0; i < n; i++) {
        vector<double> dir(n, 0.0);
        dir[i] = 1.0;
        directions.push_back(dir);
    }

    cout << "Метод Пауэлла\n";
    cout << "Целевая функция: f(x) = x1^2 + 6x2^2 + x1x2 + x1\n";
    cout << "Начальная точка: ";
    printVector(x0);
    cout << ", f = " << func(x) << endl;
    cout << "Точность: " << tol << endl << endl;
    for (int iter = 0; iter < max_iter; iter++) {
        vector<double> x_start = x;
        double f_start = func(x);
        cout << "Итерация " << iter + 1 << endl;
        cout << "Начальная точка: ";
        printVector(x);
        cout << ", f = " << f_start << endl;

        // минимизация вдоль каждого направления
        for (int i = 0; i < n; i++) {
            double lambda = golden_section(func, x, directions[i]);

            // обновление точки
            for (int j = 0; j < n; j++) {
                x[j] += lambda * directions[i][j];
            }
            cout << "Шаг " << i + 1 << ": минимизация вдоль направления ";
            printVector(directions[i]);
            cout << endl;
            cout << "    lambda = " << lambda;
            cout << ", f = " << func(x);
            cout << ", x = ";
            printVector(x);
            cout << endl;
        }

        // построение сопряженного направления
        vector<double> conjugate_dir(n);
        for (int i = 0; i < n; i++) {
            conjugate_dir[i] = x[i] - x_start[i];
        }

        double dir_norm = vector_norma(conjugate_dir);

        cout << "  Сопряженное направление: ";
        printVector(conjugate_dir);
        cout << ", норма = " << dir_norm << endl;

        if (dir_norm > tol) {
            // минимизация вдоль сопряженного направления
            vector<double> conjugate_dir_norm = normal(conjugate_dir);
            double lambda = golden_section(func, x, conjugate_dir_norm, -2.0, 2.0);

            vector<double> x_new = x;
            for (int i = 0; i < n; i++) {
                x_new[i] += lambda * conjugate_dir_norm[i];
            }

            cout << "Минимизация вдоль сопряженного направления:\n";
            cout << "    lambda = " << lambda;
            cout << ", f = " << func(x_new);
            cout << ", x = ";
            printVector(x_new);
            cout << endl;

            // обновление направлений (замена первого направления)
            directions.erase(directions.begin());
            directions.push_back(conjugate_dir_norm);
            x = x_new;
        }
        // проверка критерия остановки
        double diff_norm = 0.0;
        for (int i = 0; i < n; i++) {
            diff_norm += (x[i] - x_start[i]) * (x[i] - x_start[i]);
        }
        diff_norm = sqrt(diff_norm);
        cout << "  Изменение точки: " << diff_norm << endl;
        if (diff_norm < tol) {
            cout << "\nСходимость достигнута на итерации " << iter + 1 << " !!!" << endl;
            break;
        }
        if (iter == max_iter - 1) {
            cout << "\nДостигнуто максимальное число итераций!!!" << endl;
        }
        cout << endl;
    }

    return x;
}



int main() {
    setlocale(LC_ALL, "ru");
    // начальная точка
    vector<double> x0 = { 1.5, 1.1 };
    // запуск метода
    vector<double> result = Paul(main_function, x0, 1e-6, 100);
    // вывод результатов
    cout << "\nФинальный результат:" << endl;
    cout << "Найденный минимум: ";
    printVector(result);
    cout << endl;
    cout << "Значение функции в минимуме: f(x) = " << main_function(result) << endl;
    return 0;
}