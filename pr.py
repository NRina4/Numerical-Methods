import numpy as np


# Подынтегральная функция f(x)
def f(x):
    return 0.5 * np.cos(2 * x) * np.exp(2 * x / 5) + 2.4 * np.sin(1.5 * x) * np.exp(-6 * x) + 6 * x


# Функция для нахождения кубического корня с учетом знака
def cube_root(x):
    """Возвращает кубический корень числа с учетом его знака."""
    return np.sign(x) * (np.abs(x) ** (1 / 3))


# Функция для нахождения корней кубического уравнения x^3 + px + q = 0
def solve_cubic(p, q):
    """Находит действительные корни уравнения x^3 + px + q = 0 с помощью формулы Кардано."""
    D = q ** 2 / 4 + p ** 3 / 27
    roots = []

    if D > 0:
        # Один действительный корень
        u = cube_root(-q / 2 + np.sqrt(D))
        v = cube_root(-q / 2 - np.sqrt(D))
        root1 = u + v
        roots.append(root1)
    elif D == 0:
        # Два действительных корня, один из которых повторяется
        u = cube_root(-q / 2)
        root1 = 2 * u
        root2 = -u
        roots.extend([root1, root2])
    else:
        # Три действительных корня
        r = np.sqrt(-p ** 3 / 27)
        theta = np.arccos(-q / (2 * r))
        r_cubed_root = cube_root(-r)
        root1 = 2 * r_cubed_root * np.cos(theta / 3)
        root2 = 2 * r_cubed_root * np.cos((theta + 2 * np.pi) / 3)
        root3 = 2 * r_cubed_root * np.cos((theta + 4 * np.pi) / 3)
        roots.extend([root1, root2, root3])

    return roots


# Функция для фильтрации корней по заданному интервалу [a, b]
def filter_roots(roots, a, b):
    """Возвращает корни, которые лежат в интервале [a, b]."""
    return [root for root in roots if a <= root <= b]


# Функция для выполнения 3-точечной квадратуры Гаусса
def gauss_quadrature(f, a, b):
    """Вычисляет интеграл функции f на интервале [a, b] методом Гаусса с 3-точечной квадратурой."""
    # Стандартные узлы и веса для 3-точечной квадратуры Гаусса
    gauss_nodes = [-np.sqrt(3 / 5), 0, np.sqrt(3 / 5)]
    gauss_weights = [5 / 9, 8 / 9, 5 / 9]

    # Преобразуем узлы для интервала [a, b]
    transformed_nodes = [(b - a) / 2 * node + (a + b) / 2 for node in gauss_nodes]
    jacobian = (b - a) / 2  # Якобиан для масштабирования

    # Вычисляем интеграл
    integral = sum(gauss_weights[i] * f(transformed_nodes[i]) for i in range(3)) * jacobian
    return integral


# Основная функция для вычисления интеграла
def main():
    # Задаем интервал интегрирования
    a = 1.1
    b = 2.5

    # Вычисляем значения функции на границах интервала
    f_a = f(a)
    f_b = f(b)

    # Определяем коэффициенты для уравнения x^3 + px + q = 0
    p = -(f_a + f_b)
    q = f_a * f_b

    # Находим корни уравнения и фильтруем их по интервалу
    roots = solve_cubic(p, q)
    filtered_roots = filter_roots(roots, a, b)
    print("Найденные корни кубического уравнения на интервале [a, b]:", filtered_roots)

    # Вычисляем интеграл методом Гаусса
    integral = gauss_quadrature(f, a, b)
    print("Результат интегрирования методом Гаусса:", integral)


# Запуск основной функции
main()
