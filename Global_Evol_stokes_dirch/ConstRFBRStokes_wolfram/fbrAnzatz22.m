function fbraux = fbrAnzatz22(x1,y1,x2,y2,dtt,c,muu)
%FBRANZATZ22
%    FBRAUX = FBRANZATZ22(X1,Y1,X2,Y2,DTT,C,MUU)

%    This function was generated by the Symbolic Math Toolbox version mp('8.3').
%    07-Nov-2019 04:25:57

t2 = x1.*mp('2.0');
t3 = x1.*mp('3.0');
t4 = x2.*mp('2.0');
t5 = x2.*mp('3.0');
t6 = -x2;
t9 = -y1;
t10 = -y2;
t7 = -t4;
t8 = -t5;
t11 = t6+x1;
t12 = t10+y1;
t13 = t11.^2;
t14 = t12.^2;
t15 = t13+t14;
fbraux = -c.*mp('1.0')./(c.*t15+mp('1.0')).^(mp('9.0')./mp('2.0')).*(c.^2.*(t15.*(t11+t12).*(t9+t11+y2)+dtt.*muu.*(t3+t8+t12).*(t3+t8+t9+y2).*mp('3.0')).*mp('3.0')-c.*(t14+dtt.*muu.*mp('4.0')).*mp('3.0')+c.^3.*t15.*(t13.*t14+t13.^2.*mp('2.0')-t14.^2-dtt.*muu.*(t2+t7+t12).*(t2+t7+t9+y2).*mp('3.0'))-mp('1.0'));
