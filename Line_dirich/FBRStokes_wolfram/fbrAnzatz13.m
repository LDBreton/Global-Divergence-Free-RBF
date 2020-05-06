function fbraux = fbrAnzatz13(x1,y1,x2,y2,c,muu)
%FBRANZATZ13
%    FBRAUX = FBRANZATZ13(X1,Y1,X2,Y2,C,MUU)

%    This function was generated by the Symbolic Math Toolbox version 8.3.
%    29-Nov-2019 03:14:33

t2 = -x2;
t3 = -y2;
t4 = t2+x1;
t5 = t3+y1;
t6 = t4.^2;
t7 = t5.^2;
fbraux = c.*1.0./(c.*(t6+t7)+1.0).^(5.0./2.0).*(c.*(t6-t7.*2.0)+1.0);
