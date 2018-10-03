#include <stdio.h>
#include <math.h>
int main()
{
	double Vm,B,C,R,T,P,Tc,Pc,w,a,b,d,e,m,A;
	R=8.314;
	float V[500];
	int i=0,c,yn;
	printf ("请选择您计算的方程(前四个是计算摩尔体积后三个是计算)P：1 维里方程 2 R-K方程 3 SRK方程  4 P-R方程 5 R-K 6 SRK 7P-R");
	scanf ("%d",&c);
	if (c==1)
	{
	 printf("请依次输入 B C T P值");
	 scanf("%lf%lf%lf%lf",&B,&C,&T,&P) ;
	 V[0]=R*T/P; 
	 do{
		i++;
 		V[i]=R*T/P*(1+B/V[i-1]+C/pow(V[i-1],2));
 	    printf("第%d次迭代的值是%lf\n:",i,V[i]);
 	   } while(i<=7);
	}
	if (c==2)
	{
	 printf("请依次输入 Tc Pc T P值");
	 scanf("%lf%lf%lf%lf",&Tc,&Pc,&T,&P) ;
	 V[0]=R*T/P;
	 printf("V0是%f",V[0]) ;
	 a=0.42748*pow(R,2)*pow(Tc,2.5)/Pc;
	 b=0.08664*pow(R,1)*pow(Tc,1)/Pc;
	 printf("修正因子a是%f b是%f",a,b);
	 printf("是否采用牛顿法 :1 yes 2 no") ;
	 scanf("%d",&yn);
        if(yn==2)
	       {
		   	do{
	        	i++;
 		        V[i]=R*T/P+b-a*(V[i-1]-b)/(pow(T,0.5)*P*V[i-1]*(V[i-1]+b));
                printf("第%d次迭代的值是%lf\n:",i,V[i]);
 	           } while(i<=7);
	       }
       else
 	   {
       		do{
	        	i++;
 		        d=pow(V[i-1],3)-R*T*pow(V[i-1],2)/P-(pow(b,2)+b*R*T/P-a/P/pow(T,0.5))*V[i-1]-a*b/P/pow(T,0.5);
 		        e=3*pow(V[i-1],2)-2*R*T/P*V[i-1]-(pow(b,2)+b*R*T/P-a/P/pow(T,0.5));
 		        V[i]=V[i-1]-d/e;
                printf("第%d次迭代的值是%lf\n:",i,V[i]);
 	           } while(i<=7);
       }
	}
   if (c==3)
   {
     printf("请依次输入 Tc Pc T P w值");
	 scanf("%lf%lf%lf%lf%lf",&Tc,&Pc,&T,&P,&w) ;
	 	 V[0]=R*T/P; 
   	  m=0.48+1.574*w-0.176*pow(w,2);
   	   printf("%lf",m);
      a=0.42748*pow(R,2)*pow(Tc,2)/Pc*pow(1+m*(1-pow(T/Tc,0.5)),2);
	  b=0.08664*pow(R,1)*pow(Tc,1)/Pc;
	  printf("%lf%lf",a,b);
   	  do{
	        	i++;
	             V[i]=R*T/P+b-a*(V[i-1]-b)/(P*V[i-1]*(V[i-1]+b));
                printf("第%d次迭代的值是%lf\n:",i,V[i]); 
      } while(i<=7);
   	
   	
   }
     if (c==4)
   {
     printf("请依次输入 Tc Pc T P w值");
	 scanf("%lf%lf%lf%lf%lf",&Tc,&Pc,&T,&P,&w) ;
	 	 V[0]=R*T/P; 
   	  m=0.37464+1.54226*w-0.26992*pow(w,2);
   	   printf("%lf/n",m);
      a=0.42748*pow(R,2)*pow(Tc,2)/Pc*pow(1+m*(1-pow(T/Tc,0.5)),2);
	  b=0.07780*pow(R,1)*pow(Tc,1)/Pc;
	  printf("%lf%lf/n",a,b);
   	  do{
	        	i++;
	             V[i]=R*T/P+b-a*(V[i-1]-b)/(P*V[i-1]*(V[i-1]+b)+b);
                printf("第%d次迭代的值是%lf\n:",i,V[i]); 
      } while(i<=7);
   }
      if(c==5)
   {
  printf("请依次输入 Tc Pc T P值");
	 scanf("%lf%lf%lf",&Tc,&Pc,&T) ;
	 V[0]=0.000736;
	 printf("V0是%f",V[0]) ;
	 a=0.42748*pow(R,2)*pow(Tc,2.5)/Pc;
	 b=0.08664*pow(R,1)*pow(Tc,1)/Pc;
	 printf("修正因子a是%f b是%f",a,b);
	
 		       P =R*T/(V[0]-b)-a/(pow(T,0.5)*V[0]*(V[0]+b));
                printf("第%d次迭代的值是%lf\n:",i,P);
 	          
	   
	}
      if(c==6)
   {
   	 V[0]=0.000734;
   	 printf("请依次输入 Tc Pc w  T值");
	 scanf("%lf%lf%lf%lf",&Tc,&Pc,&w,&T) ;
	  m=0.48+1.574*w-0.176*pow(w,2);
   	   printf("%lf",m);
      a=0.42748*pow(R,2)*pow(Tc,2)/Pc*pow(1+m*(1-pow(T/Tc,0.5)),2);
	  b=0.08664*pow(R,1)*pow(Tc,1)/Pc;
	  printf("%lf%lf/n",a,b);
	  P=R*T/(V[0]-b)-a/(V[0]*(V[0]+b));
	  printf("P is %lf\n",P);
   	
   } 	
 if (c==7)
   {
     printf("请依次输入 Tc Pc T  w值");
	 scanf("%lf%lf%lf%lf",&Tc,&Pc,&T,&w) ;
	 	 V[0]=0.000734; 
   	  m=0.37464+1.54226*w-0.26992*pow(w,2);
   	   printf("%lf/n",m);
      a=0.42748*pow(R,2)*pow(Tc,2)/Pc*pow(1+m*(1-pow(T/Tc,0.5)),2);
	  b=0.07780*pow(R,1)*pow(Tc,1)/Pc;
	  printf("%lf%lf/n",a,b);
   	
	        
	             P=R*T/(V[0]-b)-a/(V[0]*(V[0]+b)+b*(V[0]-b));
                printf("第%d次迭代的值是%lf\n:",i,P); 
      
   }
   
return 0;
}
