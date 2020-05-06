function [SNormalx,SNormaly] = Symbolic_normal()

syms theta y x radiuos
assume(theta,'real')
assume(y,'real')
assume(x,'real')
assume(radiuos,'real')
f1 = radiuos*((0.8 + 0.1*(sin(6.0*theta) + sin(3.0*theta)))*cos(theta));
f2 = radiuos*((0.8 + 0.1*(sin(6.0*theta) + sin(3.0*theta)))*sin(theta)); 


Normalxaux = simplify(diff(f2,theta),'steps',5);
Normalyaux = -simplify(diff(f1,theta),'steps',5);


SNormalxaux = simplify(subs(expand(Normalxaux),theta, atan2(y,x)),'steps',5);
SNormalyaux = simplify(subs(expand(Normalyaux),theta, atan2(y,x)),'steps',5);


SNormalx = simplify(SNormalxaux/ sqrt( (SNormalxaux^2 + SNormalyaux^2) ),'steps',5);
SNormaly = simplify(SNormalyaux/ sqrt( (SNormalxaux^2 + SNormalyaux^2) ),'steps',5);

end

