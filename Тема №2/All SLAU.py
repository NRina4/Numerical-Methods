import math


class Matrix:
    def __init__(self, stroki, stolbtsy, Aq=[]):
        self.stroki = stroki
        self.stolbtsy = stolbtsy
        self.list = Aq
        if (self.list == []):
            for i in range(stroki):
                arr = []
                for j in range(stolbtsy):
                    arr.append(int(input()))
                self.list.append(arr)
        elif self.list == 'E':
            self.list = []
            for i in range(stroki):
                arr = []
                for j in range(stolbtsy):
                    if i == j:
                        arr.append(1)
                    else:
                        arr.append(0)
                self.list.append(arr)
 
    
    def pechat(self):
        for i in range(self.stroki):
            q = ''
            for j in range(self.stolbtsy):
                q+= str(self.list[i][j]) + ' '
            print(q)
    
    
    def norma1(self):
        j = 0
        summ = 0
        for i in range (self.stroki):
            summ += abs(self.list[i][j])
        mx = summ
        for j in range(1, self.stolbtsy):
            summ = 0
            for i in range (self.stroki):
                summ += abs(self.list[i][j])
            if mx < summ:
                mx = summ
        return mx


    def norma2(self):
        summa = 0
        for i in range(self.stroki):
            for j in range(self.stolbtsy):
                summa += abs(self.list[i][j])**2
        summa = summa ** 0.5
        return summa
    
    
    def normainfin(self):
        i = 0
        summ = 0
        for j in range (self.stolbtsy):
            summ += abs(self.list[i][j])
        mx = summ
        for i in range(1, self.stroki):
            summ = 0
            for j in range (self.stolbtsy):
                summ += abs(self.list[i][j])
            if mx < summ:
                mx = summ
        return mx
    

    def copy(self):
        Aq = []
        for i in range(self.stroki):
                arr = []
                for j in range(self.stolbtsy):
                    arr.append(self.list[i][j])
                Aq.append(arr)
        return Matrix(self.stroki, self.stolbtsy, Aq)


    def get_minor(self):
        Arr = []
        for i in range(1, self.stroki):
            aq = []
            for j in range(1, self.stolbtsy):
                aq.append(self.list[i][j])
            Arr.append(aq)
        return Matrix(self.stroki - 1, self.stolbtsy - 1, Arr)


    def get_column(self, j):
        j -= 1
        if j < 0 or j > self.stolbtsy:
            print("Wrong number")
        else:
            arr = []
            for i in range(self.stroki):
                a = []
                a.append(self.list[i][j])
                arr.append(a)
        return Matrix(self.stroki, 1, arr)


def multiply(m, a):
    A = m.copy().list
    for i in range(m.stroki):
            for j in range(m.stolbtsy):
                A[i][j] = round(A[i][j] * a, 12)
    return Matrix(m.stroki, m.stolbtsy, A)

def msumm(m1, m2):
    if (m1.stroki != m2.stroki) or (m1.stolbtsy != m2.stolbtsy):
        print("Wrong arguments for summ")
    else:
        A = []
        for i in range(m1.stroki):
            arr = []
            for j in range(m2.stolbtsy):
                arr.append(round(m1.list[i][j] + m2.list[i][j], 12))
            A.append(arr)
        return Matrix(m1.stroki, m1.stolbtsy, A)

def mmultiply(m1, m2):
    if m1.stolbtsy != m2.stroki:
        print("Wrong arguments for multiplication")
    else:
        Arr = []
        for i in range(m1.stroki):
            mrr = []
            for j in range(m2.stolbtsy):
                summa = 0
                k = 0
                le = 0
                while k < m1.stolbtsy:
                    while le < m2.stroki:
                        summa += m1.list[i][k] * m2.list[le][j]
                        k += 1
                        le += 1
                mrr.append(summa)
            Arr.append(mrr)
        return Matrix(m1.stroki, m2.stolbtsy, Arr)

def transpon(m):       
    A = []
    i = 0
    while i < m.stolbtsy:
        arr = []
        for j in range(m.stroki):
            arr.append(m.list[j][i])
        A.append(arr)
        i += 1
    return Matrix(m.stolbtsy, m.stroki, A)



def iteration(A0, b0, epsilon):
    A = A0.copy()
    b = b0.copy()
    u = 1 / A.normainfin()
    D = multiply(b, u)
    E = Matrix(A.stroki, A.stolbtsy, "E")
    C = msumm(E, multiply(A, -1*u))
    if C.normainfin() >= 1:
        b = mmultiply(transpon(A), b)
        A = mmultiply(transpon(A), A)
        u = 1 / A.normainfin()
        D = multiply(b, u)
        E = Matrix(A.stroki, A.stolbtsy, "E")
        C = msumm(E, multiply(A, -1*u))
    x0 = D.copy()
    iteration_count = 1
    if C.normainfin() >= 1:
        #print("Stop when |Ax-b| < e")
        x1 = msumm(mmultiply(C, x0), D).copy()
        while msumm(mmultiply(A, x1), multiply(b, -1)).normainfin() >= epsilon:
            x0 = x1.copy()
            x1 = msumm(mmultiply(C, x0), D).copy()
            iteration_count += 1
    else:
        #print("Stop when norma(B)/norma(1-B) * norma(x1 - x0) < e")
        x1 = msumm(mmultiply(C, x0), D).copy()
        delta = (C.normainfin() / (1 - C.normainfin())) * ((msumm(x1, multiply(x0, -1))).normainfin())
        while abs(delta) > epsilon:
            x0 = x1.copy()
            x1 = msumm(mmultiply(C, x0), D)
            delta = (C.normainfin() / (1 - C.normainfin())) * ((msumm(x1, multiply(x0, -1))).normainfin())
            iteration_count += 1
    print("Number of iterations is", iteration_count)
    return x1

def zeidel(A0, b0, epsilon):
    A = A0.copy()
    b = b0.copy()
    C = Matrix(A.stroki, A.stolbtsy, 'E')
    D = b.copy()
    diagonal = True
    for i in range(A.stroki):
        summa = 0
        q = abs(A.list[i][i])
        for j in range(A.stolbtsy):
            if j != i:
                summa += abs(A.list[i][j])
        if summa > q:
            diagonal = False
    if diagonal is False:
        b = mmultiply(transpon(A), b)
        A = mmultiply(transpon(A), A)
    for i in range(A.stroki):
        for j in range(A.stolbtsy):
            if i == j:
                C.list[i][j] = 0
            else:
                C.list[i][j] = -1 * A.list[i][j] / A.list[i][i]
        D.list[i][0] = b.list[i][0] / A.list[i][i]   
    iteration_count = 0
    x0 = D.copy()
    x1 = x0.copy()
    while msumm(mmultiply(A, x1), multiply(b, -1)).normainfin() >= epsilon:
        x0 = x1.copy()
        x1 = D.copy()
        for i in range(A.stroki):
            for j in range(i):
                x1.list[i][0] += C.list[i][j] * x1.list[j][0]
            for j in range(i+1, A.stolbtsy):
                x1.list[i][0] += C.list[i][j] * x0.list[j][0]
        iteration_count += 1
    print("Number of iterations is", iteration_count)
    return x1

def pLU(A0, b0):
    A = A0.copy()
    b = b0.copy()
    P = Matrix(A.stroki, A.stolbtsy, 'E')
    A2 = A.copy()
    for j in range(A.stolbtsy - 1):
        el = A2.list[j][j]
        nomer = j
        for i in range(j, A.stroki):
            if abs(A2.list[i][j]) > abs(el):
                el = A2.list[i][j]
                nomer = i
        c = A2.list[j]
        A2.list[j] = A2.list[nomer]
        A2.list[nomer] = c
        c = P.list[j]
        P.list[j] = P.list[nomer]
        P.list[nomer] = c
        
        for i in range(j + 1, A.stroki):
            A2.list[i][j] /= A2.list[j][j] 
            for k in range(j + 1, A.stolbtsy):
                A2.list[i][k] -= A2.list[i][j] * A2.list[j][k] 
    L = Matrix(A.stroki, A.stolbtsy, "E")
    U = Matrix(A.stroki, A.stolbtsy, "E")
    for i in range(A.stroki):
        for j in range(A.stolbtsy):
            if i < j:
                L.list[i][j] = 0.0
                U.list[i][j] = A2.list[i][j]
            elif i > j:
                L.list[i][j] = A2.list[i][j]
                U.list[i][j] = 0.0
            else:
                L.list[i][j] = 1.0
                U.list[i][j] = A2.list[i][j]

    b2 = mmultiply(P, b)
    x = b2.copy()
    for i in range(A.stroki):
        for j in range(i):
            x.list[i][0] -= L.list[i][j] * x.list[j][0]
    for i in range(A.stroki):
        for j in range(i):
            x.list[A.stroki - i - 1][0] -= U.list[A.stroki - i - 1][A.stolbtsy - j - 1] * x.list[A.stroki - j - 1][0]
        x.list[A.stroki - i - 1][0] /= U.list[A.stroki - i - 1][A.stolbtsy - i - 1]
    return x

def QR(A0, b0):
    A = A0.copy()
    b = b0.copy()
    Q = Matrix(A.stroki, A.stolbtsy, 'E')
    R = A.copy()
    for i in range(A.stroki - 1):
        R2 = R.copy()
        for j in range(i):
            R2 = R2.get_minor()
        y = R2.get_column(1)
        z = Matrix(R2.stroki, 1, "E")
        alpha = y.norma2()
        w = msumm(y, multiply(z, -1*alpha))
        p = w.norma2()
        w = multiply(w, 1/p)
        
        E = Matrix(R2.stroki, R2.stolbtsy, "E")
        Q1 = msumm(E, multiply(mmultiply(w, transpon(w)), -2))
        Arr = []
        for j in range(A.stroki):
            Aq = []
            for k in range(A.stolbtsy):
                if j >= i and k >= i:
                    Aq.append(Q1.list[j-i][k-i])
                elif j == k:
                    Aq.append(1)
                else:
                    Aq.append(0)
            Arr.append(Aq)
        Q1 = Matrix(A.stroki, A.stolbtsy, Arr)
        R = mmultiply(Q1, R)
        Q = mmultiply(Q, Q1)
    
    x = mmultiply(transpon(Q), b)
    for i in range(A.stroki):
        for j in range(i):
            x.list[A.stroki - i - 1][0] -= R.list[A.stroki - i - 1][A.stolbtsy - j - 1] * x.list[A.stroki - j - 1][0]
        x.list[A.stroki - i - 1][0] /= R.list[A.stroki - i - 1][A.stolbtsy - i - 1]
    return x


print('Enter the nubmer of the test (from 0 to 5) or "-1" if you want to stop')
ans = input()
while ans != '-1':
    print("Enter number of correct digits")
    epsilon = 10**(0 - int(input()))
    if ans == '0':
        A0 = Matrix(3, 3, [[0, 2, 3], [1, 2, 4], [4, 5, 6]])
        b0 = Matrix(3, 1, [[13], [17], [32]])
        
        Answer = Matrix(3, 1, [[1], [2], [3]])
        print("---------------------------------- MSI ----------------------------------")
        x = iteration(A0, b0, epsilon)
        x.pechat()
        print(" ")
        print("Answer is ")
        Answer.pechat()
        print(' ')
        print("Differrence is", msumm(Answer, multiply(x, -1)).norma2())
        print("---------------------------------- Zeidel ----------------------------------")
        x = zeidel(A0, b0, epsilon)
        x.pechat()
        print(" ")
        print("Answer is ")
        Answer.pechat()
        print(' ')
        print("Differrence is", msumm(Answer, multiply(x, -1)).norma2())
        print("---------------------------------- pLU ----------------------------------")
        x = pLU(A0, b0)
        x.pechat()
        print(" ")
        print("Answer is ")
        Answer.pechat()
        print(' ')
        print("Differrence is", msumm(Answer, multiply(x, -1)).norma2())
        print("---------------------------------- QR ----------------------------------")
        x = QR(A0, b0)
        x.pechat()
        print(" ")
        print("Answer is ")
        Answer.pechat()
        print(' ')
        print("Differrence is", msumm(Answer, multiply(x, -1)).norma2())
    elif ans == '1':
        A1 = Matrix(3, 3, [[6, 1, 1], [1, 8, 1], [1, 1, 10]])
        b1 = Matrix(3, 1, [[8], [10], [12]])
        Answer = Matrix(3, 1, [[1], [1], [1]])
        print("---------------------------------- MSI ----------------------------------")
        x = iteration(A1, b1, epsilon)
        x.pechat()
        print(" ")
        print("Answer is ")
        Answer.pechat()
        print(' ')
        print("Differrence is", msumm(Answer, multiply(x, -1)).norma2())
        print("---------------------------------- Zeidel ----------------------------------")
        x = zeidel(A1, b1, epsilon)
        x.pechat()
        print(" ")
        print("Answer is ")
        Answer.pechat()
        print(' ')
        print("Differrence is", msumm(Answer, multiply(x, -1)).norma2())
        print("---------------------------------- pLU ----------------------------------")
        x = pLU(A1, b1)
        x.pechat()
        print(" ")
        print("Answer is ")
        Answer.pechat()
        print(' ')
        print("Differrence is", msumm(Answer, multiply(x, -1)).norma2())
        print("---------------------------------- QR ----------------------------------")
        x = QR(A1, b1)
        x.pechat()
        print(" ")
        print("Answer is ")
        Answer.pechat()
        print(' ')
        print("Differrence is", msumm(Answer, multiply(x, -1)).norma2())
    elif ans == '2':
        A2 = Matrix(3, 3, [[-6, 1, 1], [1, -8, 1], [1, 1, -10]])
        b2 = Matrix(3, 1, [[-8], [-10], [-12]])
        Answer = Matrix(3, 1, [[425/227], [381/227], [353/227]])
        print("---------------------------------- MSI ----------------------------------")
        x = iteration(A2, b2, epsilon)
        x.pechat()
        print(" ")
        print("Answer is ")
        Answer.pechat()
        print(' ')
        print("Differrence is", msumm(Answer, multiply(x, -1)).norma2())
        print("---------------------------------- Zeidel ----------------------------------")
        x = zeidel(A2, b2, epsilon)
        x.pechat()
        print(" ")
        print("Answer is ")
        Answer.pechat()
        print(' ')
        print("Differrence is", msumm(Answer, multiply(x, -1)).norma2())
        print("---------------------------------- pLU ----------------------------------")
        x = pLU(A2, b2)
        x.pechat()
        print(" ")
        print("Answer is ")
        Answer.pechat()
        print(' ')
        print("Differrence is", msumm(Answer, multiply(x, -1)).norma2())
        print("---------------------------------- QR ----------------------------------")
        x = QR(A2, b2)
        x.pechat()
        print(" ")
        print("Answer is ")
        Answer.pechat()
        print(' ')
        print("Differrence is", msumm(Answer, multiply(x, -1)).norma2())
    elif ans == '3':
        A3 = Matrix(3, 3, [[-6, 7, 8], [9, -8, 5], [8, 9, -10]])
        b3 = Matrix(3, 1, [[8], [10], [12]])
        Answer = Matrix(3, 1, [[722/465], [556/465], [104/93]])
        print("---------------------------------- MSI ----------------------------------")
        x = iteration(A3, b3, epsilon)
        x.pechat()
        print(" ")
        print("Answer is ")
        Answer.pechat()
        print(' ')
        print("Differrence is", msumm(Answer, multiply(x, -1)).norma2())
        print("---------------------------------- Zeidel ----------------------------------")
        x = zeidel(A3, b3, epsilon)
        x.pechat()
        print(" ")
        print("Answer is ")
        Answer.pechat()
        print(' ')
        print("Differrence is", msumm(Answer, multiply(x, -1)).norma2())
        print("---------------------------------- pLU ----------------------------------")
        x = pLU(A3, b3)
        x.pechat()
        print(" ")
        print("Answer is ")
        Answer.pechat()
        print(' ')
        print("Differrence is", msumm(Answer, multiply(x, -1)).norma2())
        print("---------------------------------- QR ----------------------------------")
        x = QR(A3, b3)
        x.pechat()
        print(" ")
        print("Answer is ")
        Answer.pechat()
        print(' ')
        print("Differrence is", msumm(Answer, multiply(x, -1)).norma2())
    elif ans == '4':
        A4 = Matrix(3, 3, [[6, 5, 5], [5, 8, 5], [5, 5, 10]])
        b4 = Matrix(3, 1, [[8], [10], [12]])
        Answer = Matrix(3, 1, [[1/13], [9/13], [53/65]])
        print("---------------------------------- MSI ----------------------------------")
        x = iteration(A4, b4, epsilon)
        x.pechat()
        print(" ")
        print("Answer is ")
        Answer.pechat()
        print(' ')
        print("Differrence is", msumm(Answer, multiply(x, -1)).norma2())
        print("---------------------------------- Zeidel ----------------------------------")
        x = zeidel(A4, b4, epsilon)
        x.pechat()
        print(" ")
        print("Answer is ")
        Answer.pechat()
        print(' ')
        print("Differrence is", msumm(Answer, multiply(x, -1)).norma2())
        print("---------------------------------- pLU ----------------------------------")
        x = pLU(A4, b4)
        x.pechat()
        print(" ")
        print("Answer is ")
        Answer.pechat()
        print(' ')
        print("Differrence is", msumm(Answer, multiply(x, -1)).norma2())
        print("---------------------------------- QR ----------------------------------")
        x = QR(A4, b4)
        x.pechat()
        print(" ")
        print("Answer is ")
        Answer.pechat()
        print(' ')
        print("Differrence is", msumm(Answer, multiply(x, -1)).norma2())
    elif ans == '5':
        print("Insert dim")
        dim5 = int(input())
        Arr = []
        for i in range(dim5):
            Aq = []
            for j in range(dim5):
                if i > j:
                    a = 0 + 4 * epsilon
                elif i < j:
                    a = -1 - 4 * epsilon
                else:
                    a = 1 + 4 * epsilon
                Aq.append(a)
            Arr.append(Aq)    
        A5 = Matrix(dim5, dim5, Arr)
        Arr = []
        for i in range(dim5):
            Aq = []
            if i != dim5 - 1:
                Aq.append(-1)
            else:
                Aq.append(1)
            Arr.append(Aq)
        b5 = Matrix(dim5, 1, Arr)
        Arr = []
        for i in range(dim5):
            Aq = []
            if i == dim5 - 1:
                a = 1 / (1 + 4 * epsilon)
            else: 
                a = 0
            Aq.append(a)
            Arr.append(Aq)    
        Answer = Matrix(dim5, 1, Arr)
        
        print("---------------------------------- MSI ----------------------------------")
        x = iteration(A5, b5, epsilon)
        x.pechat()
        print(" ")
        print("Dim =", dim5)
        print("Answer is ")
        Answer.pechat()
        print(' ')
        print("Differrence is", msumm(Answer, multiply(x, -1)).norma2())
        print("---------------------------------- Zeidel ----------------------------------")
        x = zeidel(A5, b5, epsilon)
        x.pechat()
        print(" ")
        print("Dim =", dim5)
        print("Answer is ")
        Answer.pechat()
        print(' ')
        print("Differrence is", msumm(Answer, multiply(x, -1)).norma2())
        print("---------------------------------- pLU ----------------------------------")
        x = pLU(A5, b5)
        x.pechat()
        print(" ")
        print("Dim =", dim5)
        print("Answer is ")
        Answer.pechat()
        print(' ')
        print("Differrence is", msumm(Answer, multiply(x, -1)).norma2())
        print("---------------------------------- QR ----------------------------------")
        x = QR(A5, b5)
        x.pechat()
        print(" ")
        print("Dim =", dim5)
        print("Answer is ")
        Answer.pechat()
        print(' ')
        print("Differrence is", msumm(Answer, multiply(x, -1)).norma2())
    else:
        print("Wrong number")
        break
    print('Enter the nubmer of the test (from 0 to 5) or "-1" if you want to stop')
    ans = input()
    
