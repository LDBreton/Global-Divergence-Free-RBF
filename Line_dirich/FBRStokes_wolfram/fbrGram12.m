function fbraux = fbrGram12(x1,y1,x2,y2,c,muu)
%FBRGRAM12
%    FBRAUX = FBRGRAM12(X1,Y1,X2,Y2,C,MUU)

%    This function was generated by the Symbolic Math Toolbox version 8.3.
%    29-Nov-2019 03:13:49

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
fbraux = t2.*t5.*t6.*1.0./t11.^(1.3e+1./2.0).*(t11.^4-muu.^2.*t2.*(t10.*(t10-1.6e+1)+1.6e+1).*1.05e+2).*-3.0;