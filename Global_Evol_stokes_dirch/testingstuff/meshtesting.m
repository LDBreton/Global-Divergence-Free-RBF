addpath('FreeFem_meshing/');
Apoint = 1/10;
[P_ni,P_f] = Mesh_gen(Apoint,Apoint,1.0);
% radiuos =0.5;
% theta = 0:(2*pi/length(P_ni)):2*pi;
% x = radiuos*((0.8 + 0.1*(sin(6.0*theta) + sin(3.0*theta))).*cos(theta));
% y = radiuos*((0.8 + 0.1*(sin(6.0*theta) + sin(3.0*theta))).*sin(theta));
% P_f = [x',y'];
hold on
plot(P_ni(:,1),P_ni(:,2),'.r')
plot(P_f(:,1),P_f(:,2),'.b')
hold off
fprintf('%d %d \n',length(P_ni),length(P_f));

 
 