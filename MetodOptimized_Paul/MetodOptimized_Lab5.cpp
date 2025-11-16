#include <iostream>
#include <vector>
#include <functional>
#include <cmath>
using namespace std;

// подключение файлов с методами
#include "Paull.h"
#include "Pollak.h"

int main() {
    setlocale(LC_ALL, "ru");
    vector<double> x0 = { 1.5, 1.1 }; // начальная точка
    
    // запуск метода Паула
    vector<double> result1 = Paul(main_function, x0, 1e-6, 100);
    // вывод результатов
    cout << "\nФинальный результат:\n";
    cout << "Найденный минимум: ";
    print_vector_paul(result1);
    cout << endl;
    cout << "Значение функции в минимуме: f(x) = " << main_function(result1) << endl;
    

    // запуск метода Поллака-Райбера
    vector<double> result2 = Pollaka(x0, 1e-6, 100);
    // вывод результатов
    cout << "\nФинальные результат:\n";
    cout << "Найденный минимум: ";
    print_vector_pollak(result2);
    cout << endl;
    cout << "Значение функции в минимуме: f(x) = " << main_funct(result2) << endl;
    
    return 0;
}