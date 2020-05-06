function fbraux = fbrGram22(x1,y1,x2,y2,NX1,NY1,NX2,NY2,dtt,c,muu)
%FBRGRAM22
%    FBRAUX = FBRGRAM22(X1,Y1,X2,Y2,NX1,NY1,NX2,NY2,DTT,C,MUU)

%    This function was generated by the Symbolic Math Toolbox version 8.3.
%    05-Nov-2019 00:31:26

t2 = c.^2;
t3 = c.^3;
t5 = c.^5;
t7 = x1.*2.0;
t8 = x2.*2.0;
t9 = y1.*2.0;
t10 = y2.*2.0;
t11 = -x2;
t13 = -y2;
t4 = t2.^2;
t6 = t2.^3;
t12 = -t8;
t14 = -t10;
t15 = t11+x1;
t16 = t13+y1;
t17 = t7+t12;
t18 = t15.^2;
t19 = t9+t14;
t20 = t16.^2;
t21 = t17.^2;
t23 = t19.^2;
t25 = t18+t20;
t22 = t21.^2;
t24 = t23.^2;
t26 = c.*t25;
t27 = t26+1.0;
t28 = 1.0./t27.^(3.0./2.0);
t29 = 1.0./t27.^(5.0./2.0);
t30 = 1.0./t27.^(7.0./2.0);
t32 = 1.0./t27.^(1.1e+1./2.0);
t33 = 1.0./t27.^(1.3e+1./2.0);
t31 = t28.^3;
fbraux = c.*t28+dtt.*(c.*dtt.*t28-dtt.*t2.*t23.*t29.*(3.0./4.0))+dtt.*muu.*(t2.*t29.*1.2e+1-t3.*t21.*t30.*(1.05e+2./4.0)-t3.*t23.*t30.*(1.5e+1./4.0)+t4.*t22.*t31.*(1.05e+2./1.6e+1)+dtt.*muu.*t3.*t30.*3.6e+2+t4.*t21.*t23.*t31.*(1.05e+2./1.6e+1)-dtt.*muu.*t4.*t21.*t31.*1.575e+3-dtt.*muu.*t4.*t23.*t31.*3.15e+2+dtt.*muu.*t5.*t22.*t32.*1.0040625e+3+dtt.*muu.*t5.*t24.*t32.*(9.45e+2./1.6e+1)-dtt.*muu.*t6.*t21.^3.*t33.*1.62421875e+2+dtt.*muu.*t5.*t21.*t23.*t32.*1.063125e+3-dtt.*muu.*t6.*t21.*t24.*t33.*1.62421875e+2-dtt.*muu.*t6.*t22.*t23.*t33.*3.2484375e+2)-t2.*t21.*t29.*(3.0./4.0)+dtt.*muu.*t2.*t29.*1.2e+1-dtt.*muu.*t3.*t21.*t30.*(1.05e+2./4.0)-dtt.*muu.*t3.*t23.*t30.*(1.5e+1./4.0)+dtt.*muu.*t4.*t22.*t31.*(1.05e+2./1.6e+1)+dtt.*muu.*t4.*t21.*t23.*t31.*(1.05e+2./1.6e+1);
