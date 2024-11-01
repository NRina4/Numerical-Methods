import os
import math
import shutil

import numpy as np
import sympy as sp

import matplotlib.pyplot as plt


def my_func(x):
    """ Пользовательская функция. """
    return x ** 2 * math.cos(2 * x) + 1


def max_derivative(n, x_test):
    x = sp.Symbol('x')
    f = x ** 2 * sp.cos(2 * x) + 1
    # Вычисляем n-ую производную символически
    f_derivative = sp.diff(f, x, n)
    # Преобразуем символическую производную в числовую функцию для numpy
    f_derivative_func = sp.lambdify(x, f_derivative, 'numpy')
    # Вычисляем максимальный модуль значения производной
    max_modulus = np.abs(f_derivative_func(x_test))
    return max_modulus, f_derivative


def compute_error(method, a, b, x_nodes, x_test):
    """ Вычисление методической погрешности. """
    n = len(x_nodes)
    error = []

    if method == "evenly":
        for x in x_test:
            M, _ = max_derivative(n + 1, x)
            product = np.prod(np.abs(x - x_nodes))
            error.append(M * product / math.factorial(n + 1))
        return np.array(error)

    elif method == "chebyshev":
        for x in x_test:
            M, _ = max_derivative(n + 1, x)
            error.append((M * (b - a) ** (n + 1)) / (2 ** (2 * n + 1)) / math.factorial(n + 1))
        return np.array(error)


# Распределения ////////////////////////////////////////////////////////////////////////////////////////////////////// #
def chebyshev_distributed_nodes(a, b, n):
    """ Генерация узлов Чебышёва на интервале [a, b]. """
    i = np.arange(0, n)  # здесь n кол-во узлов, а не интервалов
    nodes = 0.5 * (b - a) * np.cos((2 * i + 1) * np.pi / (2 * n)) + 0.5 * (b + a)
    return np.sort(nodes)


def evenly_distributed_nodes(a, b, n):
    """ Генерация равномерно распределенных узлов. """
    return np.linspace(a, b, n)


def distribution(a, b, n, dist):
    """ Общая функция для распределения с выбором метода. """
    match dist:
        case "evenly":
            return evenly_distributed_nodes(a, b, n)
        case "chebyshev":
            return chebyshev_distributed_nodes(a, b, n)
        case _:
            raise ValueError("Распределение должно быть 'evenly' или 'chebyshev'")


# Интерполяции /////////////////////////////////////////////////////////////////////////////////////////////////////// #
def divided_difference(x, y, i, j):
    """ Рекурсивная функция для вычисления разделённых разностей. """
    # Используем формулу f(x_i, x_j)=(f(x_j)-f(x_i))/(x_j-x_i)
    if j == i:
        return y[i]
    elif j == i + 1:
        return (y[j] - y[i]) / (x[j] - x[i])
    else:
        return (divided_difference(x, y, i + 1, j) - divided_difference(x, y, i, j - 1)) / (x[j] - x[i])


def newton_interpolation(x, y, x_test):
    """ Интерполяция по методу Ньютона. """
    n = len(y)
    # Вычисляем разделённые разности и сохраняем их в массив
    divided_diff = [divided_difference(x, y, i, j) for i in range(n) for j in range(i, n)]

    def newton_polynomial(x_val):
        result = divided_diff[0]
        w = 1
        for i in range(1, n):
            w *= (x_val - x[i - 1])
            result += divided_diff[i] * w
        return result

    return np.array([newton_polynomial(xi) for xi in x_test])


def quadratic_spline(x_nodes, y_nodes, x_test):
    """ Квадратичный сплайн S{2,0} """
    n = len(x_nodes) - 1

    # Массивы для хранения коэффициентов многочленов
    a = np.zeros(n)
    b = np.zeros(n)
    c = np.zeros(n)

    # Находим коэффициенты a, b и c для остальных интервалов
    for i in range(0, n):
        x0, x1 = x_nodes[i], x_nodes[i + 1]
        y0, y1 = y_nodes[i], y_nodes[i + 1]

        A = np.array([
            [x0 ** 2, x0, 1],
            [x1 ** 2, x1, 1],
            [2 * x0, 1, 0]
        ])

        b_vals = np.array([y0, y1, 0])
        coeffs = np.linalg.solve(A, b_vals)
        a[i], b[i], c[i] = coeffs

    # Векторизуем вычисление значений
    y_spline = np.zeros_like(x_test, dtype=float)

    # Находим значения y_spline для каждой точки x_test
    for i in range(n):
        mask = (x_test >= x_nodes[i]) & (x_test <= x_nodes[i + 1])
        y_spline[mask] = a[i] * x_test[mask] ** 2 + b[i] * x_test[mask] + c[i]

    return y_spline


def interpolate(x_nodes, y_nodes, x_fine, method):
    """ Общая функция для интерполяции с выбором метода. """
    match method:
        case "newton":
            return newton_interpolation(x_nodes, y_nodes, x_fine)
        case "quadratic_spline":
            return quadratic_spline(x_nodes, y_nodes, x_fine)
        case _:
            raise ValueError("Метод интерполяции должен быть 'newton' или 'quadratic_spline'")


# Основные расчёты /////////////////////////////////////////////////////////////////////////////////////////////////// #
def plot_and_save_interpolations(a, b, n_values, m_values, plot_dir, interpolation_method, distribution_method):
    """ Построение графиков для каждого n с выбранным методом интерполяции. """
    for i in range(len(n_values)):
        # Генерация узлов интерполирования (n)
        x_nodes = distribution(a, b, n_values[i], distribution_method)
        x_nodes[0] = a  # необходимое условие для сплайнов (для остальных методов можно убрать, но лучше оставить)
        x_nodes[-1] = b  # необходимое условие для сплайнов (для остальных методов можно убрать, но лучше оставить)
        y_nodes = np.array([my_func(x) for x in x_nodes])

        # Генерация точек для вычисления отклонения (m)
        x_test = distribution(a, b, m_values[i], "evenly")
        y_test = np.array([my_func(xi) for xi in x_test])

        # Интерполяция
        y_interp = interpolate(x_nodes, y_nodes, x_test, interpolation_method)
        max_dev = np.max(np.abs(y_test - y_interp))

        # Вычисление методической погрешности
        error = np.max(compute_error(distribution_method, a, b, x_nodes, x_test))

        # Вывод результатов
        print(f"n = {n_values[i]}, m = {m_values[i]} ({distribution_method}): "
              f"Макс. откл. ({interpolation_method}) = {max_dev:.10f}, "
              f"Методическая погрешность = {error:.10f}")

        # Построение графика
        plt.figure(figsize=(12, 6))
        plt.plot(x_nodes, y_nodes, 'ro', label='Узлы')
        plt.plot(x_test, y_test, label='Исходная функция my_func(x)')
        plt.plot(x_test, y_interp, linestyle='--',
                 label=f'Полином {interpolation_method}, n={n_values[i]}, m={m_values[i]}')

        # Настройки графика
        plt.legend()
        plt.title(f"Интерполяция методом {interpolation_method.capitalize()}, n = {n_values[i]}, m = {m_values[i]}")
        plt.xlabel("x")
        plt.ylabel("my_func(x)")
        plt.grid(True)

        # Сохранение графика как изображения
        filename = f"{plot_dir}/plot_{distribution_method}_{interpolation_method}_{n_values[i]}_{m_values[i]}.png"
        plt.savefig(filename)
        plt.close()


# Функции для демонстрации результатов /////////////////////////////////////////////////////////////////////////////// #
def create_directory(path):
    """ Создание директории для хранения графиков. """
    if os.path.exists(path):
        shutil.rmtree(path)
    os.makedirs(path)


# //////////////////////////////////////////////////////////////////////////////////////////////////////////////////// #
def main():
    a, b = -1, 1  # Интервал [a, b]
    n_values = [5, 6, 7, 8]  # Количество узлов интерполирования
    m_values = [50, 50, 50, 50]  # Количество точек для вычисления отклонения
    print(f"Исследуемый интервал: [{a}, {b}]")

    distribution_methods = ["evenly", "chebyshev"]
    interpolation_methods = ["newton", "quadratic_spline"]

    # Путь для сохранения результатов
    plot_dir = '../plots/task4_dop'
    create_directory(plot_dir)
    for dist_method in distribution_methods:
        for inter_method in interpolation_methods:
            print()  # отступ
            plot_and_save_interpolations(a, b, n_values, m_values, plot_dir,
                                         interpolation_method=inter_method,
                                         distribution_method=dist_method)


if __name__ == "__main__":
    main()
