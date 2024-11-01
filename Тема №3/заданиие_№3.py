# f(x) = ctg x - x**2 = 0

# cos(x + 0.5) + y = 0.8
# sin y - 2x = 1.6
# __1__ #
import math
e = 10**(-4)

def f(x_):
    return math.cos(x_) / math.sin(x_) - x_**2

def fp(x_):
    return - 1 / (math.sin(x_))**2 - 2 * x_


x = math.pi - 0.001
fx0 = f(x)
while fx0 < 0:
    x -= 0.1
    fx0 = f(x)

a = x - 0.1
b = x
print(a, '\n', b, '\n')
x = a
x1 = x - f(x) / fp(x)
while abs(x1 - x) >= e:
    if a < x < b:
        x = x1
        x1 = x - f(x) / fp(x)
    else:
        x1 = (a + b) / 2
    c = f(x1)
    if c > 0:
        b = x1
    elif c < 0:
        a = x1
    else:
        print(x1)
        break
    x += e
print(x1, '\n')
print(f(x1), '\n')


# __2__ #
def f1(x__1, y__1):
    return math.cos(x__1 + 0.5) + y__1 - 0.8

def f2(x__2, y__2):
    return math.sin(y__2) - 2 * x__2 - 1.6

def f1px(x__):
    return - math.sin(x__ + 1 / 2)

def f2py(y__):
    return math.cos(y__)


f2px = - 2
f1py = 1

x = - 0.8
y = - 0.1
g = (f2(x, y) * f1py - f1(x, y) * f2py(y)) / (f1px(x) * f2py(y) - f2px * f1py)
h = (f2(x, y) * f1px(x) - f1(x, y) * f2px) / (f1py * f2px - f2py(y) * f1px(x))
x_1 = x + g
y_1 = y + h
while math.sqrt((x_1 - x)**2 + (y_1 - y)**2) >= e:
    x = x_1
    y = y_1
    g = (f2(x, y) * f1py - f1(x, y) * f2py(y)) / (f1px(x) * f2py(y) - f2px * f1py)
    h = (f2(x, y) * f1px(x) - f1(x, y) * f2px) / (f1py * f2px - f2py(y) * f1px(x))
    x_1 += g
    y_1 += h
print(x, '\n', y, '\n')

x0, y0 = 0, 0
lymbda = 0
x = (lymbda * math.sin(y0) - 1.6) / 2
y = 0.8 - lymbda * math.cos(x0 + 0.5)
g = (f2(x, y) * f1py - f1(x, y) * f2py(y)) / (f1px(x) * f2py(y) - f2px * f1py)
h = (f2(x, y) * f1px(x) - f1(x, y) * f2px) / (f1py * f2px - f2py(y) * f1px(x))
x_1 = x + g
y_1 = y + h
N = 3
for i in range(1, N + 1):
    lymbda = i / N
    x = (lymbda * math.sin(y0) - 1.6) / 2
    y = 0.8 - lymbda * math.cos(x0 + 0.5)
    g = (f2(x, y) * f1py - f1(x, y) * f2py(y) * lymbda) / (f1px(x) * f2py(y) * lymbda**2 - f2px * f1py)
    h = (f2(x, y) * f1px(x) * lymbda - f1(x, y) * f2px) / (f1py * f2px - f2py(y) * f1px(x) * lymbda**2)
    x += g
    y += h
while math.sqrt((x_1 - x)**2 + (y_1 - y)**2) >= e:
    x = x_1
    y = y_1
    g = (f2(x, y) * f1py - f1(x, y) * f2py(y)) / (f1px(x) * f2py(y) - f2px * f1py)
    h = (f2(x, y) * f1px(x) - f1(x, y) * f2px) / (f1py * f2px - f2py(y) * f1px(x))
    x_1 += g
    y_1 += h
print(x, '\n', y, '\n')
