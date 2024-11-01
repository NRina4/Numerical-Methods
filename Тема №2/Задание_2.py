import math

class Matrix:
    def __init__(self, dim, arr):
        self.dim = dim
        self.arr = arr
        if self.arr == []:
            for i in range(dim):
                m = []
                for j in range(dim):
                    m.append(int(input()))
                self.arr.append(m)
        elif self.arr == 'E':
            self.arr = []
            for i in range(dim):
                m = []
                for j in range(dim):
                    if i == j:
                        m.append(1)
                    else:
                        m.append(0)
                self.arr.append(m)
        elif self.arr == '0':
            self.arr = []
            for i in range(dim):
                m = []
                for j in range(dim):
                    m.append(0)
                self.arr.append(m)

    def print(self):
        for i in range(self.dim):
            p = ''
            for j in range(self.dim):
                p += str(self.arr[i][j]) + ' '
            print(p)

    def copy(self):
        M = []
        for i in range(self.dim):
            m = []
            for j in range(self.dim):
                m.append(self.arr[i][j])
            M.append(m)
        return Matrix(self.dim, M)

    def minor(self):
        M = []
        for i in range(1, self.dim):
            m = []
            for j in range(1, self.dim):
                m.append(self.arr[i][j])
            M. append(m)
        return Matrix(self.dim - 1, M)

    def stolb1(self):
        if self.dim >= 1:
            Y = []
            for i in range(self.dim):
                y = []
                y.append(self.arr[i][0])
                Y.append(y)
            return Vector(self.dim, Y)
        else:
            print('ошибка')

    def norm1(self):
        S = 0
        for i in range(self.dim):
            S += abs(self.arr[i][0])
        n1 = S
        for j in range(1, self.dim):
            S = 0
            for i in range(self.dim):
                S += abs(self.arr[i][j])
            n1 = S if n1 < S else n1
        return n1

    def norm2(self):
        n1 = 0
        for i in range(self.dim):
            for j in range(self.dim):
                n1 += abs(self.arr[i][j])**2
        n1 = n1**(1/2)
        return n1

    def norminf(self):
        S = 0
        for j in range(self.dim):
            S += abs(self.arr[0][j])
        ni = S
        for i in range(self.dim):
            S = 0
            for j in range(self.dim):
                S += abs(self.arr[i][j])
            ni = S if ni < S else ni
        return ni

class Vector:
    def __init__(self, dim, arr):
        self.dim = dim
        self.arr = arr
        if self.arr == []:
            for i in range(dim):
                m = []
                m.append(int(input()))
                self.arr.append(m)
        elif self.arr == '0':
            self.arr = []
            for i in range(dim):
                m = []
                m.append(0)
                self.arr.append(m)
        elif self.arr == 'E':
            self.arr = []
            for i in range(dim):
                m = []
                for j in range(1):
                    if i == 0:
                        m.append(1)
                    else:
                        m.append(0)
                self.arr.append(m)

    def print(self):
        for i in range(self.dim):
            p = ''
            p = str(self.arr[i][0])
            print(p)

    def copy(self):
        V = []
        for i in range(self.dim):
            v = []
            v.append(self.arr[i][0])
            V.append(v)
        return Vector(self.dim, V)

    def norm1(self):
        n1 = abs(self.arr[0][0])
        for j in range(1):
            S = 0
            for i in range(self.dim):
                S += abs(self.arr[i][j])
            n1 = S if n1 < S else n1
        return n1

    def norm2(self):
        S = 0
        for i in range(self.dim):
            S += (abs(self.arr[i][0]))**2
        n2 = S**0.5
        return n2

    def norminf(self):
        ni = abs(self.arr[0][0])
        for i in range(self.dim):
            S = abs(self.arr[i][0])
            ni = S if ni < S else ni
        return ni

def addition(dim, m1, m2):
    Res = []
    for i in range(dim):
        res = []
        for j in range(dim):
            res.append(m1.arr[i][j] + m2.arr[i][j])
        Res.append(res)
    return Matrix(dim, Res)

def multiplication(dim, m1, m2):
    Res = []
    for i in range(dim):
        res = []
        for j in range(dim):
            r = 0
            for k in range(dim):
                r += m1.arr[i][k] * m2.arr[k][j]
            res.append(r)
        Res.append(res)
    return Matrix(dim, Res)

def multiplication_n(dim, N, m1):
    Res = []
    for i in range(dim):
        res = []
        for j in range(dim):
            res.append(m1.arr[i][j] * N)
        Res.append(res)
    return Matrix(dim, Res)

def transp(dim, m1):
    Res = []
    for i in range(dim):
        res = []
        for j in range(dim):
            r = m1.arr[j][i]
            res.append(r)
        Res.append(res)
    return Matrix(dim, Res)

def multiplication_v(dim, v, m):
    Res = []
    for i in range(dim):
        res = []
        r = 0
        for j in range(dim):
            r += v.arr[j][0] * m.arr[i][j]
        res.append(r)
        Res.append(res)
    return Vector(dim, Res)

def addition_v(dim, v1, v2):
    Res = []
    for i in range(dim):
        res = []
        r = v1.arr[i][0] + v2.arr[i][0]
        res.append(r)
        Res.append(res)
    return Vector(dim, Res)

def multiplication_v_vt(dim, v):
    R = Matrix(dim, '0')
    for i in range(dim):
        for j in range(dim):
            R.arr[i][j] = v.arr[i][0] * v.arr[j][0]
    return R

def multiplication_v_n(dim, N, v):
    Res = []
    for i in range(dim):
        res = []
        r = 0
        for j in range(1):
            r = v.arr[i][j] * N
        res.append(r)
        Res.append(res)
    return Vector(dim, Res)


def MPI(A, b, e):
    nu = 1 / A.norminf()
    c = multiplication_v_n(b.dim, nu, b)
    B = addition(A.dim, Matrix(A.dim, 'E'), multiplication_n(A.dim, (-1) * nu, A))
    if B.norminf() >= 1:
        b = multiplication_v(b.dim, b, transp(A.dim, A))
        A = multiplication(A.dim, transp(A.dim, A), A)
        nu = 1 / A.norminf()
        c = multiplication_v_n(b.dim, nu, b)
        B = addition(A.dim, Matrix(A.dim, 'E'), multiplication_n(A.dim, (-1) * nu, A))
    k = 0
    x0 = c
    if B.norminf() >= 1:
        print('Условие остановки: |Ax - b| < e')
        x1 = addition_v(x0.dim, multiplication_v(x0.dim, x0, B), c)
        while addition_v(b.dim, multiplication_v(x1.dim, x1, A), multiplication_v_n(b.dim, -1, b)).norminf() >= e:
            x0 = x1
            x1 = addition_v(c.dim, multiplication_v(x0.dim, x0, B), c)
            k += 1
    else:
        print('Условие остановки: ||B||/(1 - ||B||)*||x^k-x^(k-1)|| < e')
        x1 = addition_v(x0.dim, multiplication_v(x0.dim, x0, B), c)
        d = (B.norminf() / (1 - B.norminf())) * ((addition_v(x1.dim, x1, multiplication_v_n(x0.dim, -1, x0))).norminf())
        while abs(d) > e:
            x0 = x1
            x1 = addition_v(x0.dim, multiplication_v(x0.dim, x0, B), c)
            d = (B.norminf() / (1-B.norminf()))*((addition_v(x1.dim, x1, multiplication_v_n(x1.dim, -1, x0))).norminf())
            k += 1
    print(k, ' - номер итерации')
    return x1

def method_Z(A, b, e):
    d = b.copy()
    C = Matrix(A.dim, 'E')
    f = 1
    for i in range(A.dim):
        s = 0
        a = A.arr[i][i]
        for j in range(A.dim):
            if i != j:
                s += abs(A.arr[i][j])
        if abs(a) < s:
            f = 0

    if f == 0:
        b = multiplication_v(b.dim, b, transp(A.dim, A))
        A = multiplication(A.dim, transp(A.dim, A), A)
    for i in range(A.dim):
        for j in range(A.dim):
            if i != j:
                C.arr[i][j] = - A.arr[i][j] / A.arr[i][i]
            else:
                C.arr[i][j] = 0
        d.arr[i][0] = b.arr[i][0] / A.arr[i][i]

    print('Условие остановки: |Ax - b| < e')
    k = 0
    x0 = d.copy()
    x1 = x0.copy()
    while (addition_v(b.dim, multiplication_v(x1.dim, x1, A), multiplication_v_n(b.dim, -1, b))).norminf() >= e:
        x0 = x1.copy()
        x1 = d.copy()
        for i in range(A.dim):
            for j in range(i):
                x1.arr[i][0] += C.arr[i][j] * x1.arr[j][0]
            for j in range(i + 1, A.dim):
                x1.arr[i][0] += C.arr[i][j] * x0.arr[j][0]
        k += 1
    print(k, ' - номер итерации')
    return x1

def LU(A, b):
    P = Matrix(A.dim, 'E')
    A2 = A.copy()
    for j in range(A2.dim - 1):
        el = A2.arr[j][j]
        eli = j
        for i in range(j, A.dim):
            if abs(A2.arr[i][j]) > abs(el):
                el = A2.arr[i][j]
                eli = i
        n = A2.arr[eli]
        A2.arr[eli] = A2.arr[j]
        A2.arr[j] = n
        np = P.arr[eli]
        P.arr[eli] = P.arr[j]
        P.arr[j] = np
        for i in range(j + 1, A2.dim):
            A2.arr[i][j] = A2.arr[i][j] / A2.arr[j][j]
            for k in range(j + 1, A2.dim):
                A2.arr[i][k] = A2.arr[i][k] - A2.arr[i][j] * A2.arr[j][k]

    L = Matrix(len(A.arr), '0')
    U = Matrix(len(A2.arr), '0')
    for i in range(A2.dim):
        for j in range(A2.dim):
            if i == j:
                L.arr[i][j] = 1
            if i <= j:
                U.arr[i][j] = A2.arr[i][j]
            else:
                L.arr[i][j] = A2.arr[i][j]

    b2 = multiplication_v(P.dim, b, P)
    y = Vector(b2.dim, '0')
    for i in range(L.dim):
        y.arr[i][0] = b2.arr[i][0]
        for j in range(i):
            y.arr[i][0] -= L.arr[i][j] * y.arr[j][0]
    x = Vector(b2.dim, '0')
    for i in range(U.dim):
        x.arr[U.dim - i - 1][0] = y.arr[U.dim - i - 1][0]
        for j in range(i):
            x.arr[U.dim - i - 1][0] -= U.arr[U.dim - i - 1][U.dim - j - 1] * x.arr[U.dim - j - 1][0]
        x.arr[U.dim - i - 1][0] /= U.arr[U.dim - i - 1][U.dim - i - 1]
    return x

def QR(A, b):
    bb = b.copy()
    Q = Matrix(A.dim, 'E')
    R = A.copy()
    for i in range(A.dim - 1):
        R0 = R.copy()
        for j in range(i):
            R0 = R0.minor()
        y = R0.stolb1()
        z = Vector(R0.dim, 'E')
        u = y.norm2()
        l = addition_v(y.dim, y, multiplication_v_n(y.dim, - u, z))
        p = l.norm2()
        w = multiplication_v_n(l.dim, 1 / p, l)

        Q0 = Matrix(R0.dim, 'E')
        Q1 = addition(Q0.dim, Q0, multiplication_n(Q0.dim, - 2, multiplication_v_vt(Q0.dim, w)))
        M = []
        for j in range(A.dim):
            m = []
            for k in range(A.dim):
                if j >= i and k >= i:
                    m.append(Q1.arr[j - i][k - i])
                elif j == k:
                    m.append(1)
                else:
                    m.append(0)
            M.append(m)
        Q1 = Matrix(A.dim, M)
        R = multiplication(Q1.dim, Q1, R)
        Q = multiplication(Q.dim, Q, Q1)

    x = multiplication_v(Q.dim, bb, transp(Q.dim, Q))
    for i in range(A.dim):
        for j in range(i):
            x.arr[A.dim - i - 1][0] -= R.arr[A.dim - i - 1][A.dim - j - 1] * x.arr[A.dim - j - 1][0]
        x.arr[A.dim - i - 1][0] /= R.arr[A.dim - i - 1][A.dim - i - 1]
    return x


print('Введите номер теста 0-5 или "выход", чтобы закончить')
n = input()
e = 10**(-8)
while n != 'выход':
    if n == '0':
        m0 = Matrix(3, [[0, 2, 3], [1, 2, 4], [4, 5, 6]])
        v0 = Vector(3, [[13], [17], [32]])

        print('Точное решение: x = (1 2 3)^(-1)')

        print('\n', '____________________MPI____________________')
        MPI(m0, v0, e).print()
        print('Абсолютная погрешность: ', addition_v(3, Vector(3, [[1], [2], [3]]), multiplication_v_n(3, -1, MPI(m0, v0, e))).norm2())

        print('\n', '____________________method_Z____________________')
        method_Z(m0, v0, e).print()
        print('Абсолютная погрешность: ', addition_v(3, Vector(3, [[1], [2], [3]]), multiplication_v_n(3, -1, method_Z(m0, v0, e))).norm2())

        print('\n', '____________________LU____________________')
        LU(m0, v0).print()
        print('Абсолютная погрешность: ', addition_v(3, Vector(3, [[1], [2], [3]]), multiplication_v_n(3, -1, LU(m0, v0))).norm2())

        print('\n', '____________________QR____________________')
        QR(m0, v0).print()
        print('Абсолютная погрешность: ', addition_v(3, Vector(3, [[1], [2], [3]]), multiplication_v_n(3, -1, QR(m0, v0))).norm2())
    elif n == '1':
        m1 = Matrix(3, [[13, 1, 1], [1, 15, 1], [1, 1, 17]])
        v1 = Vector(3, [[15], [17], [19]])

        print('Точное решение: x = (1 1 1)^(-1)')

        print('\n', '____________________MPI____________________')
        MPI(m1, v1, e).print()
        print('Абсолютная погрешность: ', addition_v(3, Vector(3, [[1], [1], [1]]), multiplication_v_n(3, -1, MPI(m1, v1, e))).norm2())

        print('\n', '____________________method_Z____________________')
        method_Z(m1, v1, e).print()
        print('Абсолютная погрешность: ', addition_v(3, Vector(3, [[1], [1], [1]]), multiplication_v_n(3, -1, method_Z(m1, v1, e))).norm2())

        print('\n', '____________________LU____________________')
        LU(m1, v1).print()
        print('Абсолютная погрешность: ', addition_v(3, Vector(3, [[1], [1], [1]]), multiplication_v_n(3, -1, LU(m1, v1))).norm2())

        print('\n', '____________________QR____________________')
        QR(m1, v1).print()
        print('Абсолютная погрешность: ', addition_v(3, Vector(3, [[1], [1], [1]]), multiplication_v_n(3, -1, QR(m1, v1))).norm2())
    elif n == '2':
        m2 = Matrix(3, [[-13, 1, 1], [1, -15, 1], [1, 1, -17]])
        v2 = Vector(3, [[-15], [-17], [-19]])

        print('Точное решение: x = (', 1105/817,  1069/817,  1041/817, ')^(-1)')

        print('/n', '____________________MPI____________________')
        MPI(m2, v2, e).print()
        print('Абсолютная погрешность: ', addition_v(3, Vector(3, [[1105/817], [1069/817], [1041/817]]), multiplication_v_n(3, -1, MPI(m2, v2, e))).norm2())

        print('/n', '____________________method_Z____________________')
        method_Z(m2, v2, e).print()
        print('Абсолютная погрешность: ', addition_v(3, Vector(3, [[1105/817], [1069/817], [1041/817]]), multiplication_v_n(3, -1, MPI(m2, v2, e))).norm2())

        print('/n', '____________________LU____________________')
        LU(m2, v2).print()
        print('Абсолютная погрешность: ', addition_v(3, Vector(3, [[1105/817], [1069/817], [1041/817]]), multiplication_v_n(3, -1, MPI(m2, v2, e))).norm2())

        print('/n', '____________________QR____________________')
        QR(m2, v2).print()
        print('Абсолютная погрешность: ', addition_v(3, Vector(3, [[1105/817], [1069/817], [1041/817]]), multiplication_v_n(3, -1, MPI(m2, v2, e))).norm2())
    elif n == '3':
        m3 = Matrix(3, [[-13, 14, 15], [16, -15, 12], [15, 16, -17]])
        v3 = Vector(3, [[15], [17], [19]])

        print('Точное решение: x = (', 8269/6362,  3559/3181,  6885/6362, ')^(-1)')

        print('\n', '____________________MPI____________________')
        MPI(m3, v3, e).print()
        print('Абсолютная погрешность: ', addition_v(3, Vector(3, [[8269/6362], [3559/3181], [6885/6362]]), multiplication_v_n(3, -1, MPI(m3, v3, e))).norm2())

        print('\n', '____________________method_Z____________________')
        method_Z(m3, v3, e).print()
        print('Абсолютная погрешность: ', addition_v(3, Vector(3, [[8269/6362], [3559/3181], [6885/6362]]), multiplication_v_n(3, -1, MPI(m3, v3, e))).norm2())

        print('\n', '____________________LU____________________')
        LU(m3, v3).print()
        print('Абсолютная погрешность: ', addition_v(3, Vector(3, [[8269/6362], [3559/3181], [6885/6362]]), multiplication_v_n(3, -1, MPI(m3, v3, e))).norm2())

        print('\n', '____________________QR____________________')
        QR(m3, v3).print()
        print('Абсолютная погрешность: ', addition_v(3, Vector(3, [[8269/6362], [3559/3181], [6885/6362]]), multiplication_v_n(3, -1, MPI(m3, v3, e))).norm2())
    elif n == '4':
        m4 = Matrix(3, [[13, 12, 12], [12, 15, 12], [12, 12, 17]])
        v4 = Vector(3, [[15], [17], [19]])

        print('Точное решение: x = (', -13/97,  181/291,  75/97, ')^(-1)')

        print('\n', '___________________MPI___________________')
        MPI(m4, v4, e).print()
        print('Абсолютная погрешность: ', addition_v(3, Vector(3, [[-13/97], [181/291], [75/97]]), multiplication_v_n(3, -1, MPI(m4, v4, e))).norm2())

        print('\n', '_________________method_Z_________________')
        method_Z(m4, v4, e).print()
        print('Абсолютная погрешность: ', addition_v(3, Vector(3, [[-13/97], [181/291], [75/97]]), multiplication_v_n(3, -1, MPI(m4, v4, e))).norm2())

        print('\n', '____________________LU____________________')
        LU(m4, v4).print()
        print('Абсолютная погрешность: ', addition_v(3, Vector(3, [[-13/97], [181/291], [75/97]]), multiplication_v_n(3, -1, MPI(m4, v4, e))).norm2())

        print('\n', '____________________QR____________________')
        QR(m4, v4).print()
        print('Абсолютная погрешность: ', addition_v(3, Vector(3, [[-13/97], [181/291], [75/97]]), multiplication_v_n(3, -1, MPI(m4, v4, e))).norm2())
    elif n == '5':
        print('Введите размерность матрицы ')
        dim_ = int(input())
        M = []
        for i in range(dim_):
            m = []
            for j in range(dim_):
                if i == j:
                    q = 1 + e * 11
                elif i > j:
                    q = 0 + e * 11
                else:
                    q = -1 - e * 11
                m.append(q)
            M.append(m)
        m5 = Matrix(dim_, M)

        M = []
        for i in range(dim_):
            m = []
            if i == (dim_ - 1):
                m.append(1)
            else:
                m.append(-1)
            M.append(m)
        v5 = Vector(dim_, M)

        M = []
        for i in range(dim_):
            m = []
            if i == (dim_ - 1):
                q = 1 / (1 + e * 11)
            else:
                q = 0
            m.append(q)
            M.append(m)
        res = Vector(dim_, M)

        print('Точное решение:')
        res.print()

        print('\n', '____________________MPI____________________')
        MPI(m5, v5, e).print()
        print('Абсолютная погрешность: ', addition_v(3, res, multiplication_v_n(3, -1, MPI(m5, v5, e))).norm2())

        print('\n', '____________________method_Z____________________')
        method_Z(m5, v5, e).print()
        print('Абсолютная погрешность: ', addition_v(3, res, multiplication_v_n(3, -1, method_Z(m5, v5, e))).norm2())

        print('\n', '____________________LU____________________')
        LU(m5, v5).print()
        print('Абсолютная погрешность: ', addition_v(3, res, multiplication_v_n(3, -1, LU(m5, v5))).norm2())

        print('\n', '____________________QR____________________')
        QR(m5, v5).print()
        print('Абсолютная погрешность: ', addition_v(3, res, multiplication_v_n(3, -1, QR(m5, v5))).norm2())
    else:
        print('Попробуйте ещё раз')
    print('Введите номер теста 0-5 или "выход", чтобы закончить')
    n = input()
