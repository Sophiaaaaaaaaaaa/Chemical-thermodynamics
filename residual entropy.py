import sympy
import math
from sympy import *
def residual_entropy(T,P,Tc,Pc,w):
    Tr=T/Tc
    Pr=P/Pc
    B0=0.083-0.422/pow(Tr,1.6)
    B1=0.139-0.172/pow(Tr,4.2)
    dB0_dTr=0.675/pow(Tr,2.6)
    dB1_dTr=0.722/pow(Tr,5.2)
    SR=-Pr*(dB0_dTr+w*dB1_dTr)*8.314
    return SR
'''
T=float(input("input T:"))29
print(type(T))
P=float(input("input P:"))
Tc=float(input("input Tc:"))
Pc=float(input("input Pc:"))
w=float(input("input w:"))
'''
'''
T1, P1, Tc1,Pc1,w1 = map(float, input("input T1, P1, Tc1,Pc1,w1 ").split())
T2,P2,Tc2,Pc2,w2=map(float, input("T2, P2, Tc2,Pc2,w2 ").split())
'''
'''
str1="294 100 283.1 5040 0.085 "
str2="509 1800 283.1 5040 0.085 "
str3="1.424 14.394e-3 -4.392e-6"
contanst1=str1.split()
contanst2=str2.split()
contanst3=str3.split()
T1=float(contanst1[0])
P1=float(contanst1[1])
Tc1=float(contanst1[2])
Pc1=float(contanst1[3])
w1=float(contanst1[4])
T2=float(contanst2[0])
P2=float(contanst2[1])
Tc2=float(contanst2[2])
Pc2=float(contanst2[3])
w2=float(contanst2[4])
A=float(contanst3[0])
B=float(contanst3[1])
C=float(contanst3[2])
print(A,B,C)
'''

#----------------输入数据，这是目前发现最简洁的写法----------------
file=open("data.txt")
a=[]
for line in file.readlines():
    a.append(line.split())
a.pop()
print(a)
T1,P1,Tc1,Pc1,w1,A,B,C=map(float,a[0][0:8])
T2,P2,Tc2,Pc2,w2,A,B,C=map(float,a[1][0:8])

#-----------------------------------------------------------------------
print(SR1,SR2)
T = symbols('T')
S1=integrate( 8.314*(A/T+B+C*T),(T,T1,T2))-8.314*log(P2/P1)
S=S1-SR1+SR2
print(S)
while ((math.fabs(S)-0)>0.00001):
    if (S<0):
        T2=(T1+T2)/2;

    else :
        T1=(T1+T2)/2
    print(T1,T2)
    SR2 = residual_entropy(T2, P2, Tc2, Pc2, w2)
    S1= integrate(8.314 * (A / T + B + C * T), (T, T1, T2)) - 8.314 * log(P2 / P1)
    S = S1 -SR1 + SR2
    print(S,SR1,SR2)
