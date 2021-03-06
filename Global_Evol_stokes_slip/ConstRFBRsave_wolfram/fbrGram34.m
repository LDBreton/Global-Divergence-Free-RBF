function fbraux = fbrGram34(x1,y1,x2,y2,NX1,NY1,NX2,NY2,dtt,c,muu)
%FBRGRAM34
%    FBRAUX = FBRGRAM34(X1,Y1,X2,Y2,NX1,NY1,NX2,NY2,DTT,C,MUU)

%    This function was generated by the Symbolic Math Toolbox version mp('8.3').
%    07-Nov-2019 13:45:12

t2 = -x2;
t3 = -y2;
t4 = t2+x1;
t5 = t3+y1;
t6 = t4.^2;
t7 = t5.^2;
t8 = t6.*mp('3.0');
t9 = t6.*mp('7.0');
t10 = t7.*mp('3.0');
t11 = t7.*mp('7.0');
t12 = -t10;
t13 = -t11;
t14 = t8+t13;
t15 = t9+t12;
t16 = c.*t14;
t17 = c.*t15;
t18 = t17+mp('2.0');
t19 = t16-mp('2.0');
t20 = NX2.*t5.*t18;
t21 = NY2.*t4.*t19;
t22 = -t20;
fbraux = c.^2.*muu.*mp('1.0')./(c.*(t6+t7)+mp('1.0')).^(mp('7.0')./mp('2.0')).*(-NX1.^2.*(t20-t21)+NY1.^2.*(t20-t21)+NX1.*NY1.*(NY2.*t5.*(c.*(t6.*mp('4.0')-t7)-mp('1.0'))+NX2.*t4.*(c.*(t6-t7.*mp('4.0'))+mp('1.0'))).*mp('4.0')).*mp('3.0');
