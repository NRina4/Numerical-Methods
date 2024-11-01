import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt

e = 1e-6


def f(x):
    """ Пользовательская функция. """
    return 0.5 * np.cos(2 * x) * np.exp(2 * x / 5) + 2.4 * np.sin(1.5 * x) * np.exp(-6 * x) + 6 * x


def p(x):
    """ Весовая функция. """
    a = 1.1
    b = 2.5
    alpha = 2 / 5
    beta = 0
    return ((x - a) ** (-alpha)) * ((b - x) ** (-beta))


def F(x):
    """ Подынтегральная функция. """
    return p(x) * f(x)


# Задание 1 ////////////////////////////////////////////////////////////////////////////////////////////////////////// #
# 1. Квадратурные формулы: Прямоугольники, Трапеция, Симпсон
def left_rectangle_method(f, a, b, n):
    h = (b - a) / n
    return h * sum(f(a + i * h) for i in range(n))


def mid_rectangle_method(f, a, b, n):
    h = (b - a) / n
    return h * sum(f(a + (i + 0.5) * h) for i in range(n))


def trapezoidal_method(f, a, b, n):
    h = (b - a) / n
    return h * (0.5 * (f(a) + f(b)) + sum(f(a + i * h) for i in range(1, n)))


def simpson_method(f, a, b, n):
    """
        Кол-во интервалов должно быть четным, так как метод основан на аппроксимации
        ф-ции с помощью парабол, которые интерполируют на группах из 3-х точек
    """
    if n % 2 == 1:
        n += 1
    h = (b - a) / n
    result = f(a) + f(b)
    for i in range(1, n, 2):  # для нечетных точек
        result += 4 * f(a + i * h)
    for i in range(2, n, 2):  # для четных точек
        result += 2 * f(a + i * h)
    return result * h / 3


# 1.2. 3-точечная формула Ньютона-Кот(е)са
def newton_cotes_3point(f, a, b, alpha, n):
    h = (b - a) / n
    result = 0
    for i in range(n):
        z1 = a + h * i
        z3 = a + h * (i + 1)
        z2 = (z1 + z3) / 2

        nu0 = ((z3 - a) ** (1 - alpha) - (z1 - a) ** (1 - alpha)) / (1 - alpha)
        nu1 = ((z3 - a) ** (2 - alpha) - (z1 - a) ** (2 - alpha)) / (2 - alpha) + nu0 * a
        nu2 = ((z3 - a) ** (3 - alpha) - (z1 - a) ** (3 - alpha)) / (3 - alpha) + nu1 * 2 * a - nu0 * a ** 2

        # matrix_z = np.array([
        #     [1, 1, 1],
        #     [z1, z2, z3],
        #     [z1 ** 2, z2 ** 2, z3 ** 2]
        # ])
        # A = np.linalg.solve(matrix_z, np.array([nu0, nu1, nu2]))

        A1 = (nu2 - nu1 * (z2 + z3) + nu0 * z2 * z3) / ((z1 - z2) * (z1 - z3))
        A2 = (nu2 - nu1 * (z1 + z3) + nu0 * z1 * z3) / ((z2 - z1) * (z2 - z3))
        A3 = (nu2 - nu1 * (z1 + z2) + nu0 * z1 * z2) / ((z3 - z1) * (z3 - z2))

        result += A1 * f(z1) + A2 * f(z2) + A3 * f(z3)
    return result


# 1.2. 3-точечная формула Гаусса
# Функция для нахождения корней Ax^3+Bx^2+Cx+D=0 методом Кардано
def cardano(A, B, C, D):
    # Приведение уравнения к каноническому виду y^3+2py+2q+0
    p = (3 * A * C - B ** 2) / (9 * A ** 2)
    q = (2 * B ** 3 - 9 * A * B * C + 27 * A ** 2 * D) / (54 * A ** 3)

    if p < 0 and (q ** 2 + p ** 3 <= 0):
        Dis = q ** 2 + p ** 3  # дискриминант
        roots = []

        if Dis < 0:
            r = np.sqrt(-p)
            # Вычисление угла phi и проверка его допустимости
            # cos_phi = q / r ** 3
            # if not -1 <= cos_phi <= 1:
            #     raise ValueError(f"Недопустимое значение для cos(phi): {cos_phi}")
            phi = np.arccos(q / r ** 3)

            # Вычисление корней
            y1 = -2 * r * np.cos(phi / 3)
            y2 = 2 * r * np.cos((np.pi - phi) / 3)
            y3 = 2 * r * np.cos((np.pi + phi) / 3)

            # Преобразование к x
            x1 = y1 - B / (3 * A)
            x2 = y2 - B / (3 * A)
            x3 = y3 - B / (3 * A)

            roots.extend([x1, x2, x3])
        elif Dis == 0:
            y = - np.sign(q) * np.cbrt(abs(q))
            x = y - B / (3 * A)

            roots.extend([x, x, x])
        else:
            u = np.cbrt(-q / 2 + np.sqrt(D))
            v = np.cbrt(-q / 2 - np.sqrt(D))
            y = u + v
            x = y - B / (3 * A)

            roots.extend([x, None, None])
        return [root for root in roots if root is not None]
    else:
        return [None, None, None]
        # raise ValueError("Уравнение не имеет трёх действительных корней.")


# Функция для фильтрации корней по заданному интервалу [a, b]
# def filter_roots(roots, a, b):
#     """Возвращает корни, которые лежат в интервале [a, b]."""
#     return [root for root in roots if a <= root <= b]


# Функция для выполнения 3-точечной квадратуры Гаусса
def gauss_3point(f, a, b, n):
    """Вычисляет интеграл функции f на интервале [a, b] методом Гаусса с 3-точечной квадратурой."""
    h = (b - a) / n
    result = 0

    for i in range(n):
        z1 = a + i * h
        z2 = z1 + h

        p = - (f(z1) + f(z2))
        q = f(z1) + f(z2)

        roots = cardano(1, 0, p, q)

        integral = 0
        for root in roots:
            if root is not None and z1 <= root <= b:
                integral += f(root)
        result += (integral )
    return result


# Задание 2 ////////////////////////////////////////////////////////////////////////////////////////////////////////// #


# 3. График зависимости погрешности от числа разбиений
def calculate_errors(method, true_value, a, b, max_n=100):
    errors = []
    ns = range(2, max_n, 2)
    for n in ns:
        integral_value = method(f, a, b, n)
        errors.append(abs(integral_value - true_value))
    return ns, errors


def main(a, b, alpha, beta, n, e):
    print(left_rectangle_method(f, a, b, n))
    print(mid_rectangle_method(f, a, b, n))
    print(trapezoidal_method(f, a, b, n))
    print(simpson_method(f, a, b, n))
    print(newton_cotes_3point(f, a, b, alpha, n))
    print(gauss_3point(f, a, b, n))

    # Точное значение интеграла для оценки погрешности
    # true_value, _ = integrate.quad(f, a, b)
    # true_value = 14.2731

    # # Построение графиков
    # plt.figure(figsize=(10, 6))
    # methods = {
    #     "Left Rectangle": left_rectangle_method,
    #     "Mid Rectangle": mid_rectangle_method,
    #     "Trapezoidal": trapezoidal_method,
    #     "Simpson": simpson_method,
    #     "Newton-Cotes 3-point": newton_cotes_3point
    # }
    #
    # for name, method in methods.items():
    #     ns, errors = calculate_errors(method, true_value, a, b)
    #     plt.plot(ns, errors, label=name)
    #
    # plt.yscale('log')
    # plt.xlabel("Number of intervals (n)")
    # plt.ylabel("Absolute error")
    # plt.title("Error dependence on the number of intervals")
    # plt.legend()
    # plt.grid(True)


if __name__ == "__main__":
    main(a=1.1, b=2.5, alpha=2 / 5, beta=0, n=100, e=1e-6)
