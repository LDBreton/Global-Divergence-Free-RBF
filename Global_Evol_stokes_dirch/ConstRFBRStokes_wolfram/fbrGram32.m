function fbraux = fbrGram32(x1,y1,x2,y2,dtt,c,muu)
%FBRGRAM32
%    FBRAUX = FBRGRAM32(X1,Y1,X2,Y2,DTT,C,MUU)

%    This function was generated by the Symbolic Math Toolbox version mp('8.3').
%    07-Nov-2019 04:20:06

t2 = -x2;
t3 = -y2;
t4 = t2+x1;
t5 = t3+y1;
t6 = t4.^2;
t7 = t5.^2;
t8 = t6+t7;
t9 = c.*t8;
fbraux = c.^2.*t4.*t5.*mp('1.0')./(t9+mp('1.0')).^(mp('9.0')./mp('2.0')).*(t9.*(t9+mp('2.0'))-c.*dtt.*muu.*(t9-mp('6.0')).*mp('5.0')+mp('1.0')).*mp('3.0');
