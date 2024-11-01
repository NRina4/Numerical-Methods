import math

e, x = 10**(-6), 0.1
dQ1 = 0.124225*e
dQ2 = 0.6882*e
dH = 0.99455*e
dC = 0.1212*e
Q1,Q2,H,C = 'Q1 = 1+(1+x)^0.5','Q2 = (1+x-x^2)^0.5','H = cosh(1+(1+x)^0.5)','C = cos((1+x-x^2)^0.5)'

print('z(x)=ch(1+(1+x)^0.5)*cos(1+x-x^2)^0.5','\n','x = [0.1, 0.2]','\n')
print('Погрешность для', Q1, ':', dQ1)
print('Погрешность для', Q2, ':', dQ2)
print('Погрешность для', H, ':', dH)
print('Погрешность для', C, ':', dC,'\n'*2)
print(' x  ','   Q1  ',' Q1(comp) ','     DQ1    ','     Q2  ',' Q2(comp)  ','    DQ2    ','     H   ','  H(comp)   ','  DH    ','     C  ','  C(comp)  ','   DC    ','     Z  ',' Z (comp)   ',' DZ ')

while x < 0.21:
    #sqrt1
    Q1 = 0
    x0 = 1
    xn = 0.5*(x0+(1+x)/x0)
    while abs(xn - x0) >= dQ1:
        Q1 = xn
        x0 = xn
        xn = 0.5*(x0+(1+x)/x0)
    Q1 += 1

    #sqrt2
    Q2 = 0
    x0 = 1
    xn = 0.5*(x0+(1+x-x**2)/x0)
    while abs(xn - x0) >= dQ2:
        Q2 = xn
        x0 = xn
        xn = 0.5*(x0+(1+x-x**2)/x0)

    #cosh
    H = 0
    k = 0
    z = (Q1)**(2*k)/ math.factorial(2*k)
    while abs(z) >= dH:
        H += z
        k += 1
        z = (Q1)**(2*k)/ math.factorial(2*k)

    #cos
    C = 0
    k = 0
    z = (-1)**k * (Q2)**(2*k) / math.factorial(2*k)
    while abs(z) >= dC:
        C += z
        k += 1
        z = (-1)**k * (Q2)**(2*k) / math.factorial(2*k)

    #z(x)
    Z = H*C
        
    print ('%.2f'%x, '%.6f'%Q1,'%.6f'%(1+math.sqrt(x+1)),'%.13f'%abs(Q1-(1+math.sqrt(x+1))),
           '%.6f'%Q2,'%.6f'%(math.sqrt(1+x-x**2)),'%.13f'%abs(Q2-(math.sqrt(1+x-x**2))),
           '%.6f'%H,'%.6f'%(math.cosh(1+math.sqrt(1+x))),'%.10f'%abs((H-(math.cosh(1+math.sqrt(1+x))))),
           '%.6f'%C,'%.6f'%(math.cos(math.sqrt(1+x-x**2))),'%.10f'%abs((C-(math.cos(math.sqrt(1+x-x**2))))),
           '%.6f'%Z,'%.6f'%(math.cosh(1+math.sqrt(1+x))*math.cos(math.sqrt(1+x-x**2))),'%.10f'%abs((Z-(math.cosh(1+math.sqrt(1+x))*math.cos(math.sqrt(1+x-x**2))))))
    x+=0.01

