addpath('FreeFem_meshing/');
addpath('locallib/');
SaveNormal = 'NormalF/';
mkdir(SaveNormal)
addpath(SaveNormal);

radiuos = 1.0;
[ P_sc, P_fc ] = Mesh_gen(1/70,1/30,radiuos);
length(P_sc)/length(P_fc)
% [NormalxSNormaly] = Symbolic_normal()
% Symbolic_normal(0.2,SaveNormal)
% Normalxvec = Normalx(P_fc(:,1),P_fc(:,2));
% Normalyvec = Normaly(P_fc(:,1),P_fc(:,2));
% theta = atan2(P_fc(:,2),P_fc(:,1));
% hold on
% plot(P_sc(:,1),P_sc(:,2),'.r')
% plot(P_fc(:,1),P_fc(:,2),'.b')
% quiver(P_fc(:,1),P_fc(:,2),Normalxvec,Normalyvec)
% %quiver(P_fc(:,1),P_fc(:,2),AngNormalx(theta),AngNormaly(theta))
% hold off

theta=(0:1/80:2*pi)';
Symbolic_normal_angular(SaveNormal)
addpath('NormalF/');
f1 = radiuos*((0.8 + 0.1*(sin(6.0*theta) + sin(3.0*theta))).*cos(theta));
f2 = radiuos*((0.8 + 0.1*(sin(6.0*theta) + sin(3.0*theta))).*sin(theta));
ANormalxvec = AngNormalx(theta,radiuos);
ANormalyvec = AngNormaly(theta,radiuos);

hold on
plot(P_sc(:,1),P_sc(:,2),'.r','MarkerSize',10)
plot(P_fc(:,1),P_fc(:,2),'+b','MarkerSize',10)
quiver(f1,f2,ANormalxvec,ANormalyvec)
grid on
hold off
xlabel('x')
ylabel('y')
grid on
legend('$\Omega$','$\partial \Omega$','Normal vector $(N_x,N_y)$','Interpreter','latex')
set(gca,'FontSize',20)
set(gca,'FontWeight','bold');