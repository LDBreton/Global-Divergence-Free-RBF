function fbraux = fbrGram11(x1,y1,x2,y2,NX1,NY1,NX2,NY2,dtt,c,muu)
%FBRGRAM11
%    FBRAUX = FBRGRAM11(X1,Y1,X2,Y2,NX1,NY1,NX2,NY2,DTT,C,MUU)

%    This function was generated by the Symbolic Math Toolbox version mp('8.3').
%    07-Nov-2019 13:07:14

t2 = c.^2;
t3 = muu.^2;
t4 = y1.*mp('2.0');
t5 = y1.*mp('3.0');
t6 = y2.*mp('2.0');
t7 = y2.*mp('3.0');
t8 = -x2;
t9 = -y2;
t10 = t8+x1;
t11 = t9+y1;
t12 = t10.^2;
t13 = t11.^2;
t14 = t12.*mp('2.0');
t15 = t13.*mp('5.0');
t17 = t12+t13;
t16 = -t15;
t18 = c.*t17;
t19 = t17.^2;
t20 = t18+mp('1.0');
fbraux = c.*mp('1.0')./t20.^(mp('1.3e+1')./mp('2.0')).*(dtt.^2.*(c.*(t14+t15-t2.*(t19.*(t12.*mp('4.0')+t16)+t3.*(t13.*mp('2.8e+1')+t17).*mp('9.0e+1')).*mp('2.0')+c.*(t3.*mp('1.8e+2')+t12.*t13.*mp('4.0')+t13.*t15-t12.^2).*mp('2.0')-c.^2.*t18.*(t19.*(t12.*mp('7.0')+t16)+t3.*(t12.*mp('1.1e+1')-t13.*mp('1.01e+2')).*mp('4.5e+1'))+t2.^2.*t19.*(t3.*(t12-t13.*mp('6.0')).*mp('4.5e+1')+t19.*(t13-t14)))+mp('1.0'))+t20.^4.*(c.*(t12-t13.*mp('2.0'))+mp('1.0'))+c.*dtt.*muu.*t20.^2.*(c.*(t5-t7+t10).*(t5-t7-x1+x2).*-mp('3.0')+t2.*t17.*(t4-t6+t10).*(t4-t6-x1+x2)+mp('4.0')).*mp('6.0'));
