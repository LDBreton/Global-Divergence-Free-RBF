function fbraux = fbrGram12(x1,y1,x2,y2,NX1,NY1,NX2,NY2,dtt,c,muu)
%FBRGRAM12
%    FBRAUX = FBRGRAM12(X1,Y1,X2,Y2,NX1,NY1,NX2,NY2,DTT,C,MUU)

%    This function was generated by the Symbolic Math Toolbox version 8.3.
%    07-Nov-2019 13:10:33

t2 = c.^2;
t3 = muu.^2;
t4 = -x2;
t5 = -y2;
t6 = t4+x1;
t7 = t5+y1;
t8 = t6.^2;
t9 = t7.^2;
t10 = t8+t9;
t11 = c.*t10;
t12 = t10.^2;
t13 = t11+1.0;
t14 = -t12;
fbraux = t2.*t6.*t7.*1.0./t13.^(1.3e+1./2.0).*(dtt.^2.*(t11.*3.0+t13-t2.*(t3.*2.8e+2+t14).*6.0+c.^2.*t11.*(t3.*4.2e+2+t12).*4.0+t2.^2.*t14.*(t3.*1.05e+2+t14))-t13.^4+c.*dtt.*muu.*t13.^2.*(t11-6.0).*1.0e+1).*-3.0;
