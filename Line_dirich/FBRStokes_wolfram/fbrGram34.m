function fbraux = fbrGram34(x1,y1,x2,y2,c,muu)
%FBRGRAM34
%    FBRAUX = FBRGRAM34(X1,Y1,X2,Y2,C,MUU)

%    This function was generated by the Symbolic Math Toolbox version 8.3.
%    29-Nov-2019 03:14:23

t2 = -x2;
t3 = -y2;
t4 = t2+x1;
t5 = t3+y1;
fbraux = c.^2.*t4.*t5.*1.0./(c.*(t4.^2+t5.^2)+1.0).^(5.0./2.0).*3.0;
