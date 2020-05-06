addpath('../Libs/AdvanpixMCT-4.5.2.12841') %% libreria de multipresicion
addpath('../Libs/DivFree_libreria/');
addpath('FreeFem_meshing/');
addpath('locallib/')


[P_ni,P_f] = Mesh_gen(1/10,1/10,1.0);
Nx = @(x,y) AngNormalx(atan2(y,x),1.0);
Ny = @(x,y) AngNormaly(atan2(y,x),1.0);

NNx = Nx(P_f(:,1),P_f(:,2));
NNy= Ny(P_f(:,1),P_f(:,2));

plot(P_ni(:,1),P_ni(:,2),'r.')
[B,I] = sort(atan2(P_f(:,2),P_f(:,1)));
quiver(P_f(:,1),P_f(:,2),NNx,NNy)