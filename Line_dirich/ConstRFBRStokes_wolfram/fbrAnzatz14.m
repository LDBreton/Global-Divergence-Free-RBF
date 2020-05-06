function fbraux = fbrAnzatz14(x1,y1,x2,y2,c,muu)
%FBRANZATZ14
%    FBRAUX = FBRANZATZ14(X1,Y1,X2,Y2,C,MUU)

%    This function was generated by the Symbolic Math Toolbox version mp('8.3').
%    29-Nov-2019 03:14:34

t2 = -x2;
t3 = -y2;
t4 = t2+x1;
t5 = t3+y1;
fbraux = c.^2.*t4.*t5.*mp('1.0')./(c.*(t4.^2+t5.^2)+mp('1.0')).^(mp('5.0')./mp('2.0')).*mp('3.0');
