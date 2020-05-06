function fbraux = fbrGram41(x1,y1,x2,y2,c,muu)
%FBRGRAM41
%    FBRAUX = FBRGRAM41(X1,Y1,X2,Y2,C,MUU)

%    This function was generated by the Symbolic Math Toolbox version 8.3.
%    29-Nov-2019 03:13:42

t2 = -x2;
t3 = -y2;
t4 = t2+x1;
t5 = t3+y1;
t6 = t4.^2;
t7 = t5.^2;
t8 = t6+t7;
t9 = c.*t8;
fbraux = c.^3.*muu.*t4.*t5.*1.0./(t9+1.0).^(9.0./2.0).*(t9-6.0).*-1.5e+1;