import numpy as np

# Подынтегральная функция f(x)
def f(x):
    return 0.5 * np.cos(2 * x) * np.exp(2 * x / 5) + 2.4 * np.sin(1.5 * x) * np.exp(-6 * x) + 6 * x

# Метод Гаусса для интегрирования с заданными весами и узлами на подотрезке
def gauss_integrate_on_interval(f, a, b):
    # Узлы и веса для 3-точечной квадратуры Гаусса на интервале [-1, 1]
    nodes = np.array([-np.sqrt(3/5), 0, np.sqrt(3/5)])
    weights = np.array([5/9, 8/9, 5/9])

    # Преобразуем узлы и веса для интервала [a, b]
    transformed_nodes = 0.5 * (b - a) * nodes + 0.5 * (b + a)
    transformed_weights = 0.5 * (b - a) * weights

    # Вычисляем интеграл с учетом узлов и весов
    integral = sum(w * f(x) for x, w in zip(transformed_nodes, transformed_weights))
    return integral

# Основная функция для вычисления интеграла на заданном интервале [a, b]
def integrate_with_gauss(f, a, b, n):
    # Размер подотрезка
    h = (b - a) / n
    total_integral = 0

    for i in range(n):
        sub_a = a + i * h
        sub_b = sub_a + h

        # Вычисляем интеграл на подотрезке с помощью метода Гаусса
        integral_value = gauss_integrate_on_interval(f, sub_a, sub_b)
        total_integral += integral_value

        # Отладочный вывод значений для каждого подотрезка
        print(f"Интеграл на подотрезке [{sub_a}, {sub_b}]: {integral_value}")

    return total_integral

# Основная функция
def main():
    a = 1.1
    b = 2.5
    n = 10  # Число подотрезков

    # Вычисляем интеграл
    integral = integrate_with_gauss(f, a, b, n)
    print("Результат интегрирования методом Гаусса:", integral)

# Запуск основной функции
main()
