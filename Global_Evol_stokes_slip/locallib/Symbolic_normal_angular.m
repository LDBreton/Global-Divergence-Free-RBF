function Symbolic_normal_angular(Dirsave)
addpath(Dirsave)
syms theta radiuos
assume(theta,'real')

f1 = radiuos*((0.8 + 0.1*(sin(6.0*theta) + sin(3.0*theta)))*cos(theta));
f2 = radiuos*((0.8 + 0.1*(sin(6.0*theta) + sin(3.0*theta)))*sin(theta)); 


Normalxaux = simplify(diff(f2,theta),'steps',5);
Normalyaux = -simplify(diff(f1,theta),'steps',5);


SNormalx = simplify(Normalxaux/ sqrt( (Normalxaux^2 + Normalyaux^2) ),'steps',5);
SNormaly = simplify(Normalyaux/ sqrt( (Normalxaux^2 + Normalyaux^2) ),'steps',5);




matlabFunction(SNormalx,'File',[Dirsave,'AngNormalx'],'vars',[theta,radiuos]);
matlabFunction(SNormaly,'File',[Dirsave,'AngNormaly'],'vars',[theta,radiuos]);
rehash;


end

