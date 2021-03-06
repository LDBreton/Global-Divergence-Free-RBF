function fbraux = fbrGram13(x1,y1,x2,y2,NX1,NY1,NX2,NY2,dtt,c,muu)
%FBRGRAM13
%    FBRAUX = FBRGRAM13(X1,Y1,X2,Y2,NX1,NY1,NX2,NY2,DTT,C,MUU)

%    This function was generated by the Symbolic Math Toolbox version 8.3.
%    07-Nov-2019 13:38:18

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
t13 = c.*t12;
fbraux = c.*(NX1.*(c.^2.*(t12.*(t8+t9).*(t8-y1+y2)-dtt.*muu.*(t3-t5+t8).*(t3-t5-x1+x2).*3.0).*3.0+c.*(t10+dtt.*muu.*4.0).*3.0+c.^2.*t13.*(t12.*(t10-t11.*2.0)+dtt.*muu.*(t2-t4+t8).*(t2-t4-x1+x2).*3.0)+1.0)+NY1.*c.*t8.*t9.*(t13.*(t13+2.0)-c.*dtt.*muu.*(t13-6.0).*5.0+1.0).*3.0).*1.0./(t13+1.0).^(9.0./2.0);
