function SNormaly = AngNormaly(theta,radiuos)
%ANGNORMALY
%    SNORMALY = ANGNORMALY(THETA,RADIUOS)

%    This function was generated by the Symbolic Math Toolbox version 8.3.
%    06-May-2020 09:29:40

t2 = theta.*3.0;
t3 = theta.*6.0;
SNormaly = radiuos.*1.0./sqrt(radiuos.^2.*(cos(theta.*9.0).*1.7e+1+cos(theta.*1.2e+1).*(3.5e+1./2.0)+cos(t2).*1.9e+1+cos(t3).*4.0+sin(t2).*1.6e+1+sin(t3).*1.6e+1+1.75e+2./2.0)).*(cos(theta.*2.0).*2.0+cos(theta.*4.0).*4.0+cos(theta.*5.0).*5.0+cos(theta.*7.0).*7.0-sin(theta).*1.6e+1).*(-1.0./2.0);
