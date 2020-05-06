function fbraux = fbrGram31(x1,y1,x2,y2,NX1,NY1,NX2,NY2,dtt,c,muu)
%FBRGRAM31
%    FBRAUX = FBRGRAM31(X1,Y1,X2,Y2,NX1,NY1,NX2,NY2,DTT,C,MUU)

%    This function was generated by the Symbolic Math Toolbox version mp('8.3').
%    07-Nov-2019 13:08:08

t2 = y1.*mp('2.0');
t3 = y2.*mp('2.0');
t4 = -x2;
t5 = -y2;
t6 = t4+x1;
t7 = t5+y1;
t8 = t6.^2;
t9 = t7.^2;
t10 = t8+t9;
t11 = c.*t10;
t12 = t11+mp('1.0');
t13 = t12.^2;
fbraux = (c.*mp('1.0')./t12.^(mp('9.0')./mp('2.0')).*(NX2.*t13.*(c.*(t8-t9.*mp('2.0'))+mp('1.0')).*mp('1.6e+1')+c.*dtt.*muu.*(NX2.*(c.*(t8-t9.*mp('9.0')).*mp('3.0')+c.*t11.*(t2-t3+t6).*(t2-t3-x1+x2)+mp('4.0')).*mp('1.6e+1')-NY2.*c.*t6.*t7.*(t11-mp('6.0')).*mp('8.0e+1')).*mp('3.0')+NY2.*c.*t6.*t7.*t13.*mp('4.8e+1')))./mp('1.6e+1');
