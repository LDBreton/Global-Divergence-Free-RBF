function [Fax1,Fax2] = interpaux(Puntos,u1,u2)
 xa = double(Puntos(:,1)) ;
 ya = double(Puntos(:,2)) ;
 v = double(u1);
 Fax1 = scatteredInterpolant(xa,ya,v,'nearest');
 Fax2 = scatteredInterpolant(xa,ya,double(u2),'nearest');

end