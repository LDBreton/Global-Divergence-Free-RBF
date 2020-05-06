function fbraux = fbrGram25(x1,y1,x2,y2,NX1,NY1,NX2,NY2,dtt,c,muu)
%FBRGRAM25
%    FBRAUX = FBRGRAM25(X1,Y1,X2,Y2,NX1,NY1,NX2,NY2,DTT,C,MUU)

%    This function was generated by the Symbolic Math Toolbox version 8.3.
%    05-Nov-2019 00:32:51

t2 = x1-x2;
t3 = y1-y2;
fbraux = c.^2.*dtt.*t2.*t3.*1.0./(c.*(t2.^2+t3.^2)+1.0).^(5.0./2.0).*-3.0;
