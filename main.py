import mpmath as mpm
from scipy.special import zeta
from sympy import totient, primerange
from mpmath import power, fdiv, exp, log, sqrt, pi, euler
from numpy import inf
from scipy.integrate import quad
from pprint import pprint

mpm.mp.dps = 25

def latex_float(f):
    float_str = format(float(f),"10.5e")
    if "e" in float_str:
        base, exponent = float_str.split("e")
        return r"${0}\cdot 10^{{{1}}}$".format(base, int(exponent))
    else:
        return float_str

def B():
    return 0.2614972128476427837554268386086958590516
def tau1(q):
    return B() + fdiv(1,2*log(q)) + fdiv(1,power(log(q),2))
def tau2(q):
    return fdiv(15 + fdiv(84.89,power(q,2/3)),8*pi) + fdiv(6 + fdiv(15,4*pi) + fdiv(3,4*pi*power(q,2/3)),log(q)) + fdiv(22.29 + fdiv(1,power(q,2/3)),power(log(q),2)) + fdiv(20.58,power(log(q),3))

print("Table 1:")
for Q in [3,4,5,10,15,20,25,50,100]:
    strng = "$" + str(Q) + "$ & $" + str(round(tau1(Q),5)) + "$ & $" + str(round(tau2(Q),5)) + "$ \\\\"
    print(strng)
print("\n")
    
print("Next, determine the constant c = 84.89...")

def eta(q):
    if q >= 599:
        return fdiv(15 + fdiv(3,totient(q)),8*pi) + fdiv(6 + fdiv(15,4*pi) + fdiv(3,4*pi*totient(q)),log(q))+ fdiv(22.29 + fdiv(1,totient(q)),power(log(q),2))+ fdiv(20.58,power(log(q),3))
    else:
        return fdiv(15 + fdiv(3,totient(q)),8*pi) + fdiv(6 + fdiv(15,4*pi) + fdiv(3,4*pi*totient(q)),log(q))+ fdiv(22.29 + fdiv(4.43 + fdiv(1,4*pi),totient(q)),power(log(q),2))+ fdiv(20.58 + fdiv(3,4*pi*totient(q)),power(log(q),3))
def eta_upper(q,c):
    return fdiv(15 + fdiv(c,power(q,2/3)),8*pi) + fdiv(6 + fdiv(15,4*pi) + fdiv(3,4*pi*power(q,2/3)),log(q)) + fdiv(22.29 + fdiv(1,power(q,2/3)),power(log(q),2)) + fdiv(20.58,power(log(q),3))

ell = 10000
while eta(ell) < eta_upper(ell,3):
    ell -= 1
print(ell)

found = False
c = 8.0
switch = True
Q = 3
while Q <= 10000:
    while eta(Q) > eta_upper(Q,c):
        c += 0.01
    Q += 1
print(c)
print("Check complete.")
print("\n")

print("Finally, compute the constant in Corollary 1.2 by taking the maximum of:")

def cor_upper_1(q):
    return 43.26016 + fdiv(sqrt(q)*abs(log(log(q))),totient(q)*log(q)) + fdiv(sqrt(q)*1.54515,totient(q)*(q-1)*log(q))
def cor_upper_2(q):
    return 43.26016 + fdiv(log(log(q)),power(q,1/6)*log(q)) + fdiv(1.54515,power(q,1/6)*(q-1)*log(q))

for Q in [j for j in range(3,31)]:
    print(Q, cor_upper_1(Q))
print(31, cor_upper_2(31))

"""
Output: 

Table 1:
$3$ & $1.54515$ & $43.26016$ \\
$4$ & $1.14251$ & $26.72431$ \\
$5$ & $0.95822$ & $19.94616$ \\
$10$ & $0.66726$ & $10.40159$ \\
$15$ & $0.58249$ & $7.92122$ \\
$20$ & $0.53983$ & $6.73171$ \\
$25$ & $0.51335$ & $6.01506$ \\
$50$ & $0.45465$ & $4.49411$ \\
$100$ & $0.41722$ & $3.58205$ \\

Next, determine the constant c = 84.89...
210
84.89000000000665
Check complete.

Finally, compute the constant in Corollary 1.2 by taking the maximum of:
3 43.94331036337837077852603
4 43.86730685393248055866763
5 43.55962422489922726412419
6 43.87003637281521482446366
7 43.46937680502737327167173
8 43.58416821787637948190339
9 43.48324569826875478289952
10 43.60546238955806033778795
11 43.40249981693180726462268
12 43.62634507509167728897217
13 43.38558379544962684609505
14 43.51755687064666965383663
15 43.45798825583965735174897
16 43.46264070540213337224797
17 43.36366495152690856112667
18 43.54205536094788171472614
19 43.3560362135257546222346
20 43.48007566365089530739681
21 43.4094998350414045014923
22 43.4425673188611413162777
23 43.3444941585180231302938
24 43.49590387164160565523918
25 43.35595535628230912232453
26 43.42226441304785808444064
27 43.36982773067104239422766
28 43.427012852522536777621
29 43.33265731077068791332321
30 43.51729940315228348752167
31 43.47132449493811542672684
"""
