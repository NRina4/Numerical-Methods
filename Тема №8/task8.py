import math
import copy
import numpy as np
import matplotlib.pyplot as plt
from dataclasses import dataclass


@dataclass
class Parameters:
    ksi: float
    A: float
    B: float
    eps: float
    p: float
    tol: float
    rtol: float
    pi: float

    c2: float
    a21: float
    b1: float
    b2: float


def func(x, y, params: Parameters):
    return np.array([params.A * y[1], - params.B * y[0]])


# Двухэтапная расчетная схема ЯМРК 2-го порядка (метод Хойна)
def runge_kutta2(x0, x1, y0, h, params: Parameters):
    n = math.ceil((x1 - x0) / h)  # runge_kutta3
    t = np.linspace(x0, x1, n + 1)
    y = np.zeros((n + 1, 2))
    y[0] = copy.deepcopy(y0)

    for i in range(n):
        k1 = h * func(t[i], y[i], params)
        k2 = h * func(t[i] + params.c2 * h, y[i] + params.a21 * k1, params)
        y[i + 1] = y[i] + params.b1 * k1 + params.b2 * k2

    return y


# Расчетная схема 3-го порядка трехэтапного ЯМРК (вторая)
def runge_kutta3(x0, x1, y0, h, params: Parameters):
    n = math.ceil((x1 - x0) / h)
    t = np.linspace(x0, x1, n + 1)
    y = np.zeros((n + 1, 2))
    y[0] = copy.deepcopy(y0)

    for i in range(n):
        k1 = h * func(t[i], y[i], params)
        k2 = h * func(t[i] + h / 3, y[i] + k1 / 3, params)
        k3 = h * func(t[i] + h * 2 / 3, y[i] + k2 * 2 / 3, params)
        y[i + 1] = y[i] + (k1 + 3 * k3) / 4

    return y


def first_step(x0, x1, y0, s, params: Parameters):
    y = copy.deepcopy(y0)

    # Пункт (a): вычисляем f(x0, y0)
    ff = func(x0, y, params)
    count = np.sum(ff == 0)  # Подсчет количества нулевых компонент

    # Пункт (b - c): вычисляем Δ и начальный шаг h
    def compute_step():
        # Вспомогательная функция для расчета шага h
        delta = (1 / max(abs(x0), abs(x1))) ** (s + 1) + np.linalg.norm(ff) ** (s + 1)
        return (params.rtol / delta) ** (1 / (s + 1))

    h1 = compute_step()

    # Пункты (d) - (e): если большинство компонент f(x0, y0) равно нулю, пересчитываем шаг
    if count > len(ff) / 2:  # Если больше половины значений f(x0, y0) нулевые
        x0 = x0 + h1
        y = y + h1 * ff  # Обновляем y по методу Эйлера
        ff = func(x0, y, params)
        h_new = compute_step()

        return min(h1, h_new)

    return h1


def full_error_rate(x0, x1, y0, h, s, params: Parameters):
    if s == 2:
        runge_kutta = runge_kutta2
    elif s == 3:
        runge_kutta = runge_kutta3
    else:
        raise ValueError("Метод поддерживает только порядок 2 и 3")

    res_1 = runge_kutta(x0, x1, y0, h, params)
    res_1 = res_1[-1]

    res_2 = runge_kutta(x0, x1, y0, h / 2, params)
    res_2 = res_2[-1]

    R = (np.array(res_2) - np.array(res_1)) / (pow(2, s) - 1)

    return np.linalg.norm(R)


def auto_step(x0, x1, y0, s, params: Parameters):
    if s == 2:
        runge_kutta = runge_kutta2
    elif s == 3:
        runge_kutta = runge_kutta3
    else:
        raise ValueError("Метод поддерживает только порядок 2 и 3")

    values = []  # Список последовательно высчитанных значений
    local_errors = []  # Список локальных погрешностей
    steps = []  # Список шагов

    x_k = x0
    y_k = copy.deepcopy(y0)
    h = first_step(x0, x1, y0, s, params)  # начальный шаг

    while x_k < x1:
        h = min(h, x1 - x_k)  # Если шаг "вылетает" за пределы, то возвращаем его

        y_curr_whole = runge_kutta(x_k, x_k + h, y_k, h, params)[-1]
        y_curr_half = runge_kutta(x_k, x_k + h, y_k, h / 2, params)[-1]

        local_error = (y_curr_half - y_curr_whole) / (1 - 2 ** -s)
        error_norm = np.linalg.norm(local_error)

        # Критерии выбора шага
        if error_norm > 2 ** s * params.p:
            h = h / 2

        elif params.p < error_norm < 2 ** s * params.p:
            x_k = x_k + h
            y_k = y_curr_half

            values.append(y_k)
            local_errors.append(error_norm)
            steps.append(h)
            h = h / 2

        elif params.p / 2 ** (s + 1) < error_norm < params.p:
            x_k = x_k + h
            y_k = y_curr_whole

            values.append(y_k)
            local_errors.append(error_norm)
            steps.append(h)

        else:
            x_k = x_k + h
            y_k = y_curr_whole

            values.append(y_k)
            local_errors.append(error_norm)
            steps.append(h)
            h = 2 * h

    return values, local_errors, steps
