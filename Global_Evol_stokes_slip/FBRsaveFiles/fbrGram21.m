function fbraux = fbrGram21(x1,y1,x2,y2,NX1,NY1,NX2,NY2,dtt,c,muu)
%FBRGRAM21
%    FBRAUX = FBRGRAM21(X1,Y1,X2,Y2,NX1,NY1,NX2,NY2,DTT,C,MUU)

%    This function was generated by the Symbolic Math Toolbox version 8.3.
%    05-Nov-2019 00:30:24

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
t21 = t17.^3;
t22 = t19.^3;
t23 = t18+t20;
t24 = c.*t23;
t25 = t24+1.0;
t26 = 1.0./t25.^(5.0./2.0);
t27 = 1.0./t25.^(7.0./2.0);
t28 = 1.0./t25.^(9.0./2.0);
t29 = 1.0./t25.^(1.1e+1./2.0);
t30 = 1.0./t25.^(1.3e+1./2.0);
fbraux = dtt.*muu.*(t3.*t17.*t19.*t27.*(4.5e+1./2.0)-t4.*t17.*t22.*t28.*(1.05e+2./1.6e+1)-t4.*t19.*t21.*t28.*(1.05e+2./1.6e+1)+dtt.*muu.*t4.*t17.*t19.*t28.*1.26e+3-dtt.*muu.*t5.*t17.*t22.*t29.*9.45e+2-dtt.*muu.*t5.*t19.*t21.*t29.*9.45e+2+dtt.*muu.*t6.*t21.*t22.*t30.*3.2484375e+2+dtt.*muu.*t6.*t17.*t19.^5.*t30.*1.62421875e+2+dtt.*muu.*t6.*t17.^5.*t19.*t30.*1.62421875e+2)+t2.*t17.*t19.*t26.*(3.0./4.0)-dtt.^2.*t2.*t17.*t19.*t26.*(3.0./4.0)+dtt.*muu.*t3.*t17.*t19.*t27.*(4.5e+1./2.0)-dtt.*muu.*t4.*t17.*t22.*t28.*(1.05e+2./1.6e+1)-dtt.*muu.*t4.*t19.*t21.*t28.*(1.05e+2./1.6e+1);
