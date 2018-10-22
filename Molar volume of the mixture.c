#include <stdio.h>
#include <math.h>
double Tc,Pc;
double pitzer (double Tc,double Pc, double T,double w)
{
double B0,B1,Tr,B;
double R=8.314;
       Tr=T/Tc;
B0=0.083-0.422/pow(Tr,1.6);
B1=0.139-0.172/pow(Tr,4.2);
B=R*Tc/Pc*(B0+w*B1);
	return B;
}
/*Prausnitz 规则  */
double Tcij (double Tci,double Tcj, double Kij)
{
	return(pow(Tci*Tcj,0.5)*(1-Kij));
}
double Vcij (double Vci,double Vcj)
{
	double Vc;
	Vc=pow((pow(Vci,0.33333)+pow(Vcj,0.3333))/2,3);
	return (Vc);
}
 double Zcij (double Zci,double Zcj)
 {
 	return ((Zci+Zcj)/2);
 }
  double Wcij (double Wci,double Wcj)
 {
 	return ((Wci+Wcj)/2);
 }
 double Pcij (double Zij,double Tij,double Vij)
 {
 	double R=8.314;
 	return(Zij*R*Tij/Vij);
 	
 }
 double Bmix (double y1,double y2,double B11,double B12,double B22)
 {
 	return (pow(y1,2)*B11+2*y1*y2*B12+pow(y2,2)*B22);
 }
 double acij (double Tij,double Pij)
 {
 	double a ;
 	double R=8.314;
 	a=0.42748*pow(R,2)*pow(Tij,2.5)/Pij;
 	return a;
 	
 }
 double amix (double y1,double y2,double a11,double a22,double a12)
 {
 	return (pow(y1,2)*a11+2*y1*y2*a12+pow(y2,2)*a22);
 }
  double bci (double Ti,double Pi)
 {
 	double b ;
 	double R=8.314;
 	b=0.08664*R*Ti/Pi;
 	return b;
 	
 }
  double bmix (double y1,double y2,double b1,double b2)
 {
 	return (y1*b1+y2*b2);
 }
 double pbhjc (double y1,double y2,double a,double b)
{
	return (y1*a+y2*b);
}
/*普遍化关联图的二元体系混合计算*/ 
 int main()
 {
 	double y1,y2,T,P,Tc1,Tc2,Pc1,Pc2,Vc1,Vc2,Zc1,Zc2,w1,w2,k12;
 	double Tc12,Vc12,Zc12,Pc12,w;
 	double a11,a22,a12,a;
 	double b1,b2,b;
 	double B11,B22,B12,B,B1;
 	double Zmix[20],hmix[20]; /*用于R-K方程 迭代的Z*/ 
 	double Z;/*用于普通化关联的维里系数*/ 
 	double R=8.314;
 	double Pmr,Tmr,wm;/*用于普遍化关联图，参数查找*/ 
 	double Vm1,Vm2,Vm3;
 	double H,S,H1,S1;
 	double c1,c2,c3,c4,c5,c6;
 	int i=0;
 	double fi1,fi2,derta,fia,fib;
 	double f1,f2,fa,fb;
 	double rk1,rk2;
 	/* R-K方程*/ 
 	printf("请依次输入状态函数:");
 	scanf("%lf%lf",&T,&P);
 	printf("请输入两组分各自的百分数:\n");
 	scanf("%lf%lf",&y1,&y2);
 	printf("请输入两组分各自的Tc:\n");
 	scanf("%lf%lf",&Tc1,&Tc2);
 	printf("请输入两组分各自的Pc:\n");
 	scanf("%lf%lf",&Pc1,&Pc2);
	printf("请输入两组分各自的Vc:\n");
 	scanf("%lf%lf",&Vc1,&Vc2);
	printf("请输入两组分各自的Zc:\n");
 	scanf("%lf%lf",&Zc1,&Zc2);
 	printf("请输入两组分各自的w:\n");
 	scanf("%lf%lf",&w1,&w2);
 	printf("请输入两组分之间的K12:\n");
 	scanf("%lf",&k12);
 	Tc12=Tcij(Tc1,Tc2,k12);
 	Vc12=Vcij(Vc1,Vc2);
 	Zc12=Zcij(Zc1,Zc2);
 	Pc12=Pcij(Zc12,Tc12,Vc12);
 	printf("Tc12,Vc12,Pc12 is %lf %lf %lf\n",Tc12,Vc12,Pc12);
 	w=Wcij (w1,w2);
 	a11=acij(Tc1,Pc1);
 	a12=acij(Tc2,Pc2);
 	a22=acij(Tc12,Pc12);
 	a=amix(y1,y2,a11,a12,a22);
 	printf("a11 a12 a22 amix is %lf %lf %lf %lf\n",a11,a12,a22,a) ;
 	b1=bci(Tc1,Pc1);
 	b2=bci(Tc2,Pc2);
 	b=bmix(y1,y2,b1,b2);
 	printf("bmix is %lf",b) ;
 	Zmix[0]=1;
 	do
 	{
	 	hmix[i]=b*P/(Zmix[i]*R*T);
	 	i++;
	 	Zmix[i]=1/(1-hmix[i-1])-a/(b*R*pow(T,1.5))*(hmix[i-1]/(1+hmix[i-1]));
	 }while(i>=10);
	 	printf("zmix is %lf\n",Zmix[i]) ;
	 Vm1=Zmix[i]*R*T/P;
	 rk1=log(Vm1/(Vm1-b))+b1/(Vm1-b)-2*(y1*a11+y2*a12)/b/R/pow(T,1.5)*log(Vm1+b/Vm1)+a*b1/(pow(b,2)*R*pow(T,1.5))*(log((Vm1+b)/Vm1)-(b/(Vm1+b)))-log(Zmix[i]);
	 rk2=log(Vm1/(Vm1-b))+b1/(Vm1-b)-2*(y1*a12+y2*a22)/b/R/pow(T,1.5)*log(Vm1+b/Vm1)+a*b2/(pow(b,2)*R*pow(T,1.5))*(log((Vm1+b)/Vm1)-(b/(Vm1+b)))-log(Zmix[i]);
	 fa=pow(2.71828,rk1)*P;
	 fb=pow(2.71828,rk2)*P;
	 printf("fa fb is %lf %lf\n",fa,fb);
	 
	 
 /*普遍化维里系数*/
 	
 	B11=pitzer(Tc1,Pc1,T,w1);
 	B22=pitzer(Tc2,Pc2,T,w2);
 	B12=pitzer(Tc12,Pc12,T,w) ;
 	B=Bmix(y1,y2,B11,B22,B12);
 	printf("B11 B22 B12 Bmix is%lf%lf%lf%lf\n",B11,B22,B12,B );
 	derta=2*B12-B11-B22;
 	fi1=pow(2.71828,(P/R/T*(B11+y1*y1*derta)));
 	fi2=pow(2.71828,(P/R/T*(B22+y2*y2*derta)));
 	printf("derta fi1 fi2 is %lf%lf%lf\n",derta,fi1,fi2);
 	f1=fi1*y1*P;
 	f2=fi2*y2*P;
 	printf("f1 f2 is %lf  %lf\n",f1,f2);
 	Z=1+B*P/R/T;
 	printf("Z is %lf\n",Z) ;
 	Vm2=Z*R*T/P;	
 	printf("由R-K方程计算得来的Vm 是:%lf\n",Vm1);
 	printf("由普遍化维里系数计算得来的Vm 是:%lf\n",Vm2);
 	
 	H=R*T*((-1.5)*(a/b/R/pow(T,1.5))*log(1+b*P/Z/R/T)+Z-1);
 	S=log(Z*(1-b*P/Z/R/T))-0.5*a/b/R/pow(T,1.5)*log(1+b*P/Z/R/T);
 	printf("由R-K方程计算得来的剩余焓 剩余熵：H is %lf,S is %lf:\n",H,S);
    c1=0.675/pow(T/Tc1,2.6);
    c2=0.722/pow(T/Tc1,5.2);
    c3=0.675/pow(T/Tc2,2.6);
    c4=0.722/pow(T/Tc2,5.2);
    c5=0.675/pow(T/Tc12,2.6);
    c6=0.722/pow(T/Tc12,5.2);
    B1=y1*y1*R/Pc1*(c1+w1*c2)+2*y1*y2*R/Pc12*(c5+w*c6)+y2*y2*R/Pc2*(c3+w2*c4);
    H1=P*T*(B/T-B1);
    S1=-P*B1;
    printf("由普遍化维里方程计算得来的剩余焓 剩余熵：H1 is %lf,S1 is %lf:\n" ,H1,S1);
    Pmr=P/pbhjc(y1,y2,Pc1,Pc2);
    Tmr=T/pbhjc(y1,y2,Tc1,Tc2);
    wm=pbhjc(y1,y2,w1,w2);
    printf("Pmr,Tmr,w is %lf %lf %lf",Pmr,Tmr,wm);
 return 0;	
}
