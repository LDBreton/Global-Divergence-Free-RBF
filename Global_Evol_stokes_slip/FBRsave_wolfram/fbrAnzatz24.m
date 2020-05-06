function fbraux = fbrAnzatz24(x1,y1,x2,y2,NX1,NY1,NX2,NY2,dtt,c,muu)
%FBRANZATZ24
%    FBRAUX = FBRANZATZ24(X1,Y1,X2,Y2,NX1,NY1,NX2,NY2,DTT,C,MUU)

%    This function was generated by the Symbolic Math Toolbox version 8.3.
%    07-Nov-2019 13:46:30

t2 = -x2;
t3 = -y2;
t4 = t2+x1;
t5 = t3+y1;
t6 = t4.^2;
t7 = t5.^2;
t8 = t6.*3.0;
t9 = t7.*7.0;
t10 = -t9;
t11 = t8+t10;
t12 = c.*t11;
t13 = t12-2.0;
fbraux = c.^2.*muu.*1.0./(c.*(t6+t7)+1.0).^(7.0./2.0).*(NX2.^2.*t4.*t13-NY2.^2.*t4.*t13+NX2.*NY2.*t5.*(c.*(t6.*4.0-t7)-1.0).*4.0).*-3.0;
