function SNormaly = AngNormaly(theta,radiuos)
%ANGNORMALY
%    SNORMALY = ANGNORMALY(THETA,RADIUOS)

%    This function was generated by the Symbolic Math Toolbox version mp('8.3').
%    11-Nov-2019 00:44:02

t2 = theta.*mp('3.0');
t3 = theta.*mp('6.0');
SNormaly = radiuos.*mp('1.0')./sqrt(radiuos.^2.*(cos(theta.*mp('9.0')).*mp('1.7e+1')+cos(theta.*mp('1.2e+1')).*(mp('3.5e+1')./mp('2.0'))+cos(t2).*mp('1.9e+1')+cos(t3).*mp('4.0')+sin(t2).*mp('1.6e+1')+sin(t3).*mp('1.6e+1')+mp('1.75e+2')./mp('2.0'))).*(cos(theta.*mp('2.0')).*mp('2.0')+cos(theta.*mp('4.0')).*mp('4.0')+cos(theta.*mp('5.0')).*mp('5.0')+cos(theta.*mp('7.0')).*mp('7.0')-sin(theta).*mp('1.6e+1')).*(-mp('1.0')./mp('2.0'));
