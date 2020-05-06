function fbraux = fbrGram22(x1,y1,x2,y2,c,muu)
%FBRGRAM22
%    FBRAUX = FBRGRAM22(X1,Y1,X2,Y2,C,MUU)

%    This function was generated by the Symbolic Math Toolbox version 8.3.
%    29-Nov-2019 03:14:12

t2 = c.^2;
t3 = -x2;
t4 = -y2;
t5 = t3+x1;
t6 = t4+y1;
t7 = t5.^2;
t8 = t6.^2;
t9 = t7+t8;
t10 = c.*t9;
t11 = t10+1.0;
fbraux = -c.*1.0./t11.^(1.3e+1./2.0).*(-t11.^5+c.*t8.*t11.^4.*3.0+muu.^2.*t2.*(c.*(t7.*2.8e+1+t9).*4.0-t2.*t9.*(t7.*1.01e+2-t8.*1.1e+1)+c.*t10.^2.*(t7.*6.0-t8)-8.0).*4.5e+1);