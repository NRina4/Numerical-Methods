import numpy as np


def mnt(z1, z3, a, alpha):
    """Вычисление моментов. """
    nu0 = ((z3 - a) ** (1 - alpha) - (z1 - a) ** (1 - alpha)) / (1 - alpha)
    nu1 = ((z3 - a) ** (2 - alpha) - (z1 - a) ** (2 - alpha)) / (2 - alpha) + nu0 * a
    nu2 = ((z3 - a) ** (3 - alpha) - (z1 - a) ** (3 - alpha)) / (3 - alpha) + nu1 * 2 * a - nu0 * a ** 2
    return np.array([nu0, nu1, nu2])


def mnt_new(i, z1, z3, alpha):
    """Вычисление моментов порядка i на каждом частичном промежутке [a, b] в новых переменных (t = x - a)"""
    return (z1**(i-alpha+1) - z3**(i-alpha+1)) / (i-alpha+1)


# Задание 1.1 //////////////////////////////////////////////////////////////////////////////////////////////////////// #
# Квадратурные формулы: Прямоугольники, Трапеция, Симпсон
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


# Задание 1.2 //////////////////////////////////////////////////////////////////////////////////////////////////////// #
# 3-точечная формула Ньютона-Кот(е)са
def newton_cotes_3point(f, a, b, n, alpha):
    h = (b - a) / n

    result = 0
    for i in range(n-1):
        z1 = a + i * h
        z3 = a + (i + 1) * h
        z2 = (z1 + z3) / 2

        nu = mnt(z1, z3, a, alpha)

        # Собственные подсчёты из учебника
        # A1 = (nu[2] - nu[1] * (z2 + z3) + nu[0] * z2 * z3) / ((z1 - z2) * (z1 - z3))
        # A2 = (nu[2] - nu[1] * (z1 + z3) + nu[0] * z1 * z3) / ((z2 - z1) * (z2 - z3))
        # A3 = (nu[2] - nu[1] * (z1 + z2) + nu[0] * z1 * z2) / ((z3 - z1) * (z3 - z2))
        # A = [A1, A2, A3]

        z = np.array([z1, z2, z3])
        c = np.array([np.ones(3), z, z ** 2])
        A = np.linalg.solve(c, nu)

        result += A[0] * f(z1) + A[1] * f(z2) + A[2] * f(z3)
    return result


# 1.2. 3-точечная формула Гаусса
# Функция для нахождения корней Ax^3+Bx^2+Cx+D=0 методом Кардано
def cardano(A, B, C, D):
    # Приведение уравнения к каноническому виду y^3+2py+2q+0
    p = (3 * A * C - B ** 2) / (9 * A ** 2)
    q = (2 * (B ** 3) - 9 * A * B * C + 27 * (A ** 2) * D) / (54 * A ** 3)
    Dis = q ** 2 + p ** 3

    if p < 0 and Dis <= 0:
        if q < 0:
            r = (-1) * np.sqrt(abs(p))
        else:
            r = np.sqrt(abs(p))

        # Вычисление угла phi и проверка его допустимости
        phi = np.arccos(q / r ** 3)

        # Вычисление корней
        y1 = -2 * r * np.cos(phi / 3)
        y2 = 2 * r * np.cos((np.pi - phi) / 3)
        y3 = 2 * r * np.cos((np.pi + phi) / 3)

        # Преобразование к x
        x1 = y1 - B / (3 * A)
        x2 = y2 - B / (3 * A)
        x3 = y3 - B / (3 * A)

        return np.array([x1, x2, x3])
    else:
        raise ValueError("Уравнение не имеет трёх действительных корней.")


# Функция для выполнения 3-точечной квадратуры Гаусса
def gauss_3point(f, a, b, n, alpha):
    a_old , b_old = a, b
    # Замена t = x - a
    b = b_old - a
    a = a_old - a
    h = (b - a) / n

    result = 0
    for i in range(n - 1):
        z1 = a + i * h
        z3 = a + (i + 1) * h
        z2 = (z1 + z3) / 2

        # 1) посчитать моменты
        mu0 = mnt_new(0, z1, z3, alpha)
        mu1 = mnt_new(1, z1, z3, alpha)
        mu2 = mnt_new(2, z1, z3, alpha)
        mu3 = mnt_new(3, z1, z3, alpha)
        mu4 = mnt_new(4, z1, z3, alpha)
        mu5 = mnt_new(5, z1, z3, alpha)

        # 2) вычислить A_j через СЛАУ
        matr = np.array([
            [mu0, mu1, mu2],
            [mu1, mu2, mu3],
            [mu2, mu3, mu4]
        ])
        great_moments = - np.array([mu3, mu4, mu5])
        Aj = np.linalg.solve(matr, great_moments)

        # 3) найти узлы x_j как корни кубического ур-я с коэффициентами A_j из (2) — по формулам Кардано
        xj = cardano(1, Aj[2], Aj[1], Aj[0])

        # 4) найти A
        c = np.array([np.ones(3), xj, xj ** 2])
        mu_vec = - np.array([mu0, mu1, mu2])
        A = np.linalg.solve(c, mu_vec)
        result += A[0] * f(z1 + a_old) + A[1] * f(z2 + a_old) + A[2] * f(z3 + a_old)
    return result


# Задание 2 ////////////////////////////////////////////////////////////////////////////////////////////////////////// #

