function fbraux = fbrGram43(x1,y1,x2,y2,NX1,NY1,NX2,NY2,dtt,c,muu)
%FBRGRAM43
%    FBRAUX = FBRGRAM43(X1,Y1,X2,Y2,NX1,NY1,NX2,NY2,DTT,C,MUU)

%    This function was generated by the Symbolic Math Toolbox version 8.3.
%    05-Nov-2019 00:32:37

t2 = c.^2;
t3 = c.^3;
t4 = x1.*2.0;
t5 = x2.*2.0;
t6 = y1.*2.0;
t7 = y2.*2.0;
t8 = -x2;
t10 = -y2;
t9 = -t5;
t11 = -t7;
t12 = t8+x1;
t13 = t10+y1;
t14 = t4+t9;
t15 = t12.^2;
t16 = t6+t11;
t17 = t13.^2;
t18 = t14.^2;
t19 = t14.^3;
t20 = t16.^2;
t21 = t16.^3;
t22 = t15+t17;
t23 = c.*t22;
t24 = t23+1.0;
t25 = 1.0./t24.^(5.0./2.0);
t26 = 1.0./t24.^(7.0./2.0);
t27 = t2.*t14.*t25.*3.0;
t28 = t2.*t16.*t25.*3.0;
t29 = t3.*t19.*t26.*(1.5e+1./8.0);
t30 = t3.*t21.*t26.*(1.5e+1./8.0);
t33 = t3.*t14.*t20.*t26.*(1.5e+1./4.0);
t34 = t3.*t16.*t18.*t26.*(1.5e+1./4.0);
t37 = t3.*t14.*t20.*t26.*(1.5e+1./8.0);
t38 = t3.*t16.*t18.*t26.*(1.5e+1./8.0);
t31 = -t29;
t32 = -t30;
t35 = -t33;
t36 = -t34;
t39 = t27+t35;
t40 = t28+t36;
t41 = t27+t31+t37;
t42 = t28+t32+t38;
fbraux = NX1.*(NX2.*(NX2.*muu.*t42-NY2.*muu.*t39)-NY2.*(NX2.*muu.*t39+NY2.*muu.*t42))+NY1.*(NX2.*(NX2.*muu.*t41+NY2.*muu.*t40)+NY2.*(NX2.*muu.*t40-NY2.*muu.*t41));