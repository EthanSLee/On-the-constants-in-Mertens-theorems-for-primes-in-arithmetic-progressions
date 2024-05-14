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
    return 2*B() + fdiv(1,log(q)) + fdiv(2,power(log(q),2))
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
    return 43.26016 + fdiv(sqrt(q)*abs(log(log(q))),totient(q)*log(q)) + fdiv(sqrt(q)*3.0903,totient(q)*(q-1)*log(q))
def cor_upper_2(q):
    return 43.26016 + fdiv(log(log(q)),power(q,1/6)*log(q)) + fdiv(3.0903,power(q,1/6)*(q-1)*log(q))

for Q in [j for j in range(3,31)]:
    print(Q, cor_upper_1(Q))
print(31, cor_upper_2(31))

"""
Output: 

Table 1:
$3$ & $3.0903$ & $43.26016$ \\
$4$ & $2.28503$ & $26.72431$ \\
$5$ & $1.91644$ & $19.94616$ \\
$10$ & $1.33451$ & $10.40159$ \\
$15$ & $1.16498$ & $7.92122$ \\
$20$ & $1.07966$ & $6.73171$ \\
$25$ & $1.02669$ & $6.01506$ \\
$50$ & $0.9093$ & $4.49411$ \\
$100$ & $0.83445$ & $3.58205$ \\


Next, determine the constant c = 84.89...
210
84.89000000000665
Check complete.


Finally, compute the constant in Corollary 1.2 by taking the maximum of:
3 44.55232373709540948736947
4 44.23883689433741086578603
5 43.69379608424608751875047
6 44.08127170551821699043981
7 43.52773399645205478260895
8 43.65922862101490468626737
9 43.52719745267570658977721
10 43.66440813398191746654016
11 43.42387140414394219719294
12 43.67530023962002499725155
13 43.40066725751865309062179
14 43.5456429576373614198773
15 43.47771891262811470718382
16 43.48121720742237988760389
17 43.37244861853181919152572
18 43.56429117440571753395693
19 43.36309613512320594084758
20 43.49525102558833996801925
21 43.41919041730827356245065
22 43.45373229301537560218679
23 43.34937711887736237437034
24 43.50884872678725448028336
25 43.36095564625176112338642
26 43.4303251011515653658721
27 43.37503297318144443447418
28 43.43458594052991659237654
29 43.33580920989800077959842
30 43.52802475263904119954004
31 43.47978683427913607520876
"""
