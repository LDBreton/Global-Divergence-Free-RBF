function fbraux = fbrAnzatz11(x1,y1,x2,y2,dtt,c,muu)
%FBRANZATZ11
%    FBRAUX = FBRANZATZ11(X1,Y1,X2,Y2,DTT,C,MUU)

%    This function was generated by the Symbolic Math Toolbox version 8.3.
%    07-Nov-2019 04:25:43

t2 = y1.*2.0;
t3 = y1.*3.0;
t4 = y2.*2.0;
t5 = y2.*3.0;
t6 = -x2;
t7 = -y2;
t8 = t6+x1;
t9 = t7+y1;
t10 = t8.^2;
t11 = t9.^2;
t12 = t10+t11;
fbraux = c.*1.0./(c.*t12+1.0).^(9.0./2.0).*(c.^2.*(t12.*(t8+t9).*(t8-y1+y2)-dtt.*muu.*(t3-t5+t8).*(t3-t5-x1+x2).*3.0).*3.0+c.*(t10+dtt.*muu.*4.0).*3.0+c.^3.*t12.*(t12.*(t10-t11.*2.0)+dtt.*muu.*(t2-t4+t8).*(t2-t4-x1+x2).*3.0)+1.0);
