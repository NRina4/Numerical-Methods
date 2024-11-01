import numpy as np
import sympy as sp

# Определяем переменную x как символ для использования в sympy


def max_derivative(n, x_test):
    x = sp.Symbol('x')
    f = x ** 2 * sp.cos(2 * x) + 1
    # Вычисляем n-ую производную символически
    f_derivative = sp.diff(f, x, n)
    # Преобразуем символическую производную в числовую функцию для numpy
    f_derivative_func = sp.lambdify(x, f_derivative, 'numpy')
    # Находим значения производной в заданных точках и вычисляем максимальный модуль
    max_modulus = np.max(np.abs(f_derivative_func(x_test)))
    return max_modulus, f_derivative


a, b = 0, 10
n = 4
points = np.linspace(a, b, 10)

f_derivative, max_modulus = max_derivative(n+1, points)
# Вывод результата
print(f"Производная порядка {n + 1}: {f_derivative}")
print(f"Максимум модуля производной на заданных точках: {max_modulus}")