function fbraux = fbrGram15(x1,y1,x2,y2,NX1,NY1,NX2,NY2,dtt,c,muu)
%FBRGRAM15
%    FBRAUX = FBRGRAM15(X1,Y1,X2,Y2,NX1,NY1,NX2,NY2,DTT,C,MUU)

%    This function was generated by the Symbolic Math Toolbox version 8.3.
%    07-Nov-2019 13:46:04

t2 = -x2;
t3 = -y2;
t4 = t2+x1;
t5 = t3+y1;
t6 = t4.^2;
t7 = t5.^2;
fbraux = -c.*dtt.*(c.*(t6.*2.0-t7)-1.0).*1.0./(c.*(t6+t7)+1.0).^(5.0./2.0);
