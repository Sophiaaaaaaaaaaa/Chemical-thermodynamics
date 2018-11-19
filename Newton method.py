import sympy
import math
from matplotlib import pyplot as plt
from sympy import *
x=Symbol('x')
#函数输入区
f=x*x*x
N=x-f/diff(f)
#初值设置
b=[]
b.append(2)
b.append(N.evalf(subs = {x:b[0]}))
i=1
#循环区
while  math.fabs(b[i]-b[i-1])>0.001:
    b.append(N.evalf(subs = {x:b[i]}))
    i+=1
g=[]
print(b)

for e in range(i+1):
    g.append(e)

plt.plot(g, b, color='#4B0082', linestyle='--', linewidth=2, alpha=0.5)
plt.show()
