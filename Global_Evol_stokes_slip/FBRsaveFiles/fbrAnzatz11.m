function fbraux = fbrAnzatz11(x1,y1,x2,y2,NX1,NY1,NX2,NY2,dtt,c,muu)
%FBRANZATZ11
%    FBRAUX = FBRANZATZ11(X1,Y1,X2,Y2,NX1,NY1,NX2,NY2,DTT,C,MUU)

%    This function was generated by the Symbolic Math Toolbox version 8.3.
%    05-Nov-2019 00:32:53

t2 = c.^2;
t3 = c.^3;
t5 = x1.*2.0;
t6 = x2.*2.0;
t7 = y1.*2.0;
t8 = y2.*2.0;
t9 = -x2;
t11 = -y2;
t4 = t2.^2;
t10 = -t6;
t12 = -t8;
t13 = t9+x1;
t14 = t11+y1;
t15 = t5+t10;
t16 = t13.^2;
t17 = t7+t12;
t18 = t14.^2;
t19 = t15.^2;
t20 = t17.^2;
t21 = t16+t18;
t22 = c.*t21;
t23 = t22+1.0;
t24 = 1.0./t23.^(5.0./2.0);
t25 = 1.0./t23.^(7.0./2.0);
t26 = 1.0./t23.^(9.0./2.0);
fbraux = c.*1.0./t23.^(3.0./2.0)-t2.*t20.*t24.*(3.0./4.0)+dtt.*muu.*t2.*t24.*1.2e+1-dtt.*muu.*t3.*t19.*t25.*(1.5e+1./4.0)-dtt.*muu.*t3.*t20.*t25.*(1.05e+2./4.0)+dtt.*muu.*t4.*t20.^2.*t26.*(1.05e+2./1.6e+1)+dtt.*muu.*t4.*t19.*t20.*t26.*(1.05e+2./1.6e+1);
