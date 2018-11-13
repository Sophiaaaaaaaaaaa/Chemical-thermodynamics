import math
global gamma1,gamma2,ps1,ps2,phi1s,phi2s,p,y1,y2,b1,b2,d12,y3,y4,ps4;
global tc12,vc12,zc12,w12,pc12,tr1,tr2,tr12,b_0_1,b_0_2,b_0_12,b_1_1,b_1_2,b_1_12,b1,b2,b12,d12,ps
ps=1
tc12=1
vc12=1
zc12=1
w12=1
pc12=1
tr1=1
tr2=1
tr12=1
b_0_1=1
b_0_2=1
b_0_12=1
b_1_1=1
b_1_2=1
b_1_12=1
b1=1
b2=1
b12=1
d12=1
y3=0
y4=0
ps4=0
d12=1
gamma1=1
gamma2=1
ps1=1
ps2=1
y1=1
y2=1
b1=1
b2=1
p=1000
#常数输入区
phi1s=1;
phi2s=1;
#str1="yz pz T x1 x2"
str1="0 0 311.65 0.4 0.6"
# str2= "EPS A1 B1 C1 A2 B2 C2 Pc1 Pc2 Tc1 Tc2 Vc1 Vc2 Zc1 Zc2 W1  W2 "
str2="0.01 	16.4918 	3593.39 	-35.2249 	14.1603 	2948.78 	-44.5633 	8100000 	4890000 	512.6 	562.1	0.000118 0.000259 0.224 	0.271 	0.559	0.212 	 "
#str3="M_A12 M_A21 V_A12 V_A21 W_A12 W_A21 Vml Vm2 N_A12 N_A21 "
str3="2.1873 	1.6654 	2.2350 	1.6871 	7187.6873	823.5087  0.00004073 0.00008941	3351.6647 	4813.0347 	 "
#str4="ALPHA12 q1 q2 r1 r2 U_A12 U_A21 R"
str4="0.5011 	1.4320 	2.4000 	1.4311 	3.1878	-310.3963 	4639.0050 	8.314"
constant1=str1.split();
constant2=str2.split();
constant3=str3.split();
constant4=str4.split();
#读取数据
#各个环境变量
EPS=float(constant2[0])
A1=float(constant2[1])
B1=float(constant2[2])
C1=float(constant2[3])
A2=float(constant2[4])
B2=float(constant2[5])
C2=float(constant2[6])
Pc1=float(constant2[7])
Pc2=float(constant2[8])
Tc1=float(constant2[9])
Tc2=float(constant2[10])
Vc1=float(constant2[11])
Vc2=float(constant2[12])
Zc1=float(constant2[13])
Zc2=float(constant2[14])
W1=float(constant2[15])
W2=float(constant2[16])

M_A12=float(constant3[0])
M_A21=float(constant3[1])
V_A12=float(constant3[2])
V_A21=float(constant3[3])
W_A12=float(constant3[4])
W_A21=float(constant3[5])
Vm1=float(constant3[6])
Vm2=float(constant3[7])
N_A12=float(constant3[8])
N_A21=float(constant3[9])

ALPHA12	= float(constant4[0])
q1	= float(constant4[1])
q2	= float(constant4[2])
r1	= float(constant4[3])
r2	= float(constant4[4])
U_A12	= float(constant4[5])
U_A21	= float(constant4[6])
R	= float(constant4[7])

#输入环境变量
yz	=float(constant1[0])
pz	=float(constant1[1])
T	=float(constant1[2])
x1	=float(constant1[3])
x2	=float(constant1[4])

def antoine():
  global ps1,ps2
  ps1=math.exp(A1-B1/(T+C1))*1000;
  ps2=math.exp(A2-B2/(T+C2))*1000;
  x2=1.0-x1;

def margules():
    global gamma1,gamma2
    ln_gamma1 = x2 * x2 * (M_A12 + 2 * (M_A21 - M_A12) * x1);
    ln_gamma2 = x1 * x1 * (M_A21 + 2 * (M_A12 - M_A21) * x2);
    gamma1 = math.exp(ln_gamma1);
    gamma2 = math.exp(ln_gamma2);

def  van_laar():
    global gamma1, gamma2
    ln_gamma1 = V_A12 * math.pow((1 + V_A12 / V_A21 * x1 / x2), -2)
    ln_gamma2 = V_A21 * math.pow((1 + V_A21 / V_A12 * x2 / x1), -2)
    gamma1 = math.exp(ln_gamma1);
    gamma2 = math.exp(ln_gamma2);

def  wilson():
    global gamma1, gamma2
    w_A12 = Vm2 / Vm1 * math.exp(-W_A12 / R / T)
    w_A21 = Vm1 / Vm2 * math.exp(-W_A21 / R / T)
    ln_gamma1 = -math.log(x1 + w_A12 * x2)+x2*(w_A12/(x1+w_A12*x2)-w_A21/(x2+w_A21*x1))
    ln_gamma2 = -math.log(x2 + w_A21 * x1) + x1 * (w_A21 / (x2+w_A21 * x1) - w_A12 / (x1 + w_A12 * x2))
    gamma1 = math.exp(ln_gamma1);
    gamma2 = math.exp(ln_gamma2);

def NRTL():
    global gamma1, gamma2
    t12 = N_A12 / R / T
    t21 = N_A21 / R / T;
    g12 = math.exp(-ALPHA12 * t12);
    g21 = math.exp(-ALPHA12 * t21);
    ln_gamma1 = x2 * x2 * (t21 * g21 * g21 / (x1 + x2 * g21) / (x1 + x2 * g21) + t12 * g12 / (x2 + x1 * g12) / (x2 + x1 * g12))
    ln_gamma2=x1 * x1 * (t12 * g12 * g12 / (x2+x1 * g12) / (x2 + x1 * g12) + t21 * g21 / (x1 + x2 * g21) / (x1 + x2 *g21))
    gamma1=math.exp(ln_gamma1);
    gamma2=math.exp(ln_gamma2);

def UNIQUAC():
    global gamma1, gamma2
    sita1 = q1 * x1 / (q1 * x1 + q2 * x2);
    sita2 = q2 * x2 / (q1 * x1 + q2 * x2);
    phi1 = r1 * x1 / (r1 * x1 + r2 * x2);
    phi2 = r2 * x2 / (r1 * x1 + r2 * x2)
    t12 = math.exp(-U_A12 / R / T);
    t21 = math.exp(-U_A21 / R / T);
    l1 = 5 * (r1 - q1) - (r1- 1);
    l2 = 5 * (r2 - q2) - (r2 - 1);
    ln_gamma1 = math.log(phi1 / x1) + 5 * q1 * math.log(sita1 / phi1) + phi2 * (l1 - r1 / r2 * l2) - q1 * math.log(
        sita1 + sita2 * t21) + sita2 * q1 * (t21 / (sita1 + sita2 * t21) - t12 / (sita2 + sita1 * t12))
    ln_gamma2 = math.log(phi2 / x2) + 5 * q2 * math.log(sita2 / phi2) + phi1 * (l2 - r2 / r1 * l1) - q2 * math.log(
        sita2 + sita1 * t12) + sita1 * q2 * (t12 / (sita2 + sita1 * t12) - t21 / (sita1 + sita2 * t21))

    gamma1 = math.exp(ln_gamma1);
    gamma2 = math.exp(ln_gamma2);

def virial():
    global tc12, vc12, zc12, w12, pc12, tr1, tr2, tr12, b_0_1, b_0_2, b_0_12, b_1_1, b_1_2, b_1_12, b1, b2, b12,d12
    tc12 = float(math.sqrt(Tc1 * Tc2))
    vcl2 = math.pow(((math.pow(Vc1, 1 / 3) + math.pow(Vc2, 1 / 3)) / 2), 3)
    zc12 = float((Zc1 + Zc2) / 2)
    w12 = float((W1 + W2) / 2)
    pc12 = zc12* R * tc12 / vcl2
    tr1 = T / Tc1;
    tr2 = T / Tc2;
    tr12 = T / tc12;
    b_0_1 = 0.083 - 0.422 / pow(tr1, 1.6)
    b_1_1 = 0.139 - 0.172 / pow(tr1, 4.2);
    b_0_2 = 0.083 - 0.22 / pow(tr2, 1.6);
    b_1_2 = 0.139 - 0.172 / pow(tr2, 4.2)
    b_0_12 = 0.083 - 0.422 / pow(tr12, 1.6);
    b_1_12 = 0.139 - 0.172 / pow(tr12, 4.2);
    b1 = R * Tc1 / Pc1 * (b_0_1 + W1 * b_1_1)
    b2 = R * Tc2 / Pc2 * (b_0_2 + W2 * b_1_2)
    b12 = R * tc12 / pc12 * (b_0_12 + w12 * b_1_12)
    d12 = 2 * b12 - b1 - b2;
def calc():
 global gamma1,gamma2,ps1,ps2,phi1s,phi2s,p,ps,b1,b2,d12,y1,y2,y3,y4,ps4

 ps= x1 * gamma1 * ps1 /phi1s+ x2*gamma2 * ps2/phi2s


 while math.fabs((p - ps) ) > EPS:
      p = ps;
      y1 = x1 * gamma1 * ps1 / (phi1s*p);
      y2 = x2 * gamma2 * ps2 / ( phi2s*p);
      y1 = y1 / (y1 + y2);
      y2 = 1 - y1;
      phi1s =math.exp((b1 * (p - ps1) + p * y2 * y2 * d12) / (R * T))
      phi2s = math.exp((b2 * (p - ps2) + p * y1 * y1 * d12) / (R * T))
      ps = x1 * gamma1 * ps1/phi1s  + x2 * gamma2 * ps2/phi2s ;


 y3=y1
 y4=y2
 ps4=ps


print(EPS,M_A12,M_A21,A1,A2)
antoine()

virial()

margules()
print(gamma1,gamma2)
print("Results with Margules Eqn.:");
calc()
print(y3,y4,ps4,'\n');
van_laar()
print(V_A12,V_A21,gamma1,gamma2)
print("Results with van_laar() Eqn.:");
calc()
print(y3,y4,ps4,'\n');
wilson()
print(gamma1,gamma2)
print("Results with wilson Eqn.:");
calc()
print(y3,y4,ps4,'\n');
NRTL()
print(gamma1,gamma2)
print("Results with NRTL Eqn.:");
calc()
print(y3,y4,ps4,'\n');
UNIQUAC()
print(gamma1,gamma2)
print("Results with UNIQUAC Eqn.:");
calc()
print(y3,y4,ps4,'\n');
