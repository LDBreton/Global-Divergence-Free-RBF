function fbraux = fbrGram27(x1,y1,x2,y2,dtt,c,muu)
%FBRGRAM27
%    FBRAUX = FBRGRAM27(X1,Y1,X2,Y2,DTT,C,MUU)

%    This function was generated by the Symbolic Math Toolbox version mp('8.3').
%    07-Nov-2019 04:24:12

t2 = -x2;
t3 = -y2;
t4 = t2+x1;
t5 = t3+y1;
t6 = t4.^2;
t7 = t5.^2;
t8 = t6+t7;
t9 = c.*t8;
fbraux = c.^3.*muu.*t4.*t5.*mp('1.0')./(t9+mp('1.0')).^(mp('1.3e+1')./mp('2.0')).*(-t9.*(t9.*(t9-mp('4.0'))-mp('1.1e+1'))+c.*dtt.*muu.*(t9.*(t9-mp('1.6e+1'))+mp('1.6e+1')).*mp('2.1e+1')+mp('6.0')).*mp('1.5e+1');