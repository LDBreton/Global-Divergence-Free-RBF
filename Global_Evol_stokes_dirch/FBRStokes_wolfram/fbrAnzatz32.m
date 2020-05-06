function fbraux = fbrAnzatz32(x1,y1,x2,y2,dtt,c,muu)
%FBRANZATZ32
%    FBRAUX = FBRANZATZ32(X1,Y1,X2,Y2,DTT,C,MUU)

%    This function was generated by the Symbolic Math Toolbox version 8.3.
%    07-Nov-2019 04:26:02

t2 = -y2;
t3 = t2+y1;
fbraux = c.*dtt.*t3.*1.0./(c.*((x1-x2).^2+t3.^2)+1.0).^(3.0./2.0);