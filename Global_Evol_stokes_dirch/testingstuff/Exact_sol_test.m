%Funcion para obtener:
%f1,f2,dirchi1,dirchi2,p_exact_x,p_exact_y
%Usnado la solucion exacta de Wendland

function [f1,f2,u_exact1,u_exact2,p_x,p_y,L1,L2,ut_exact1,ut_exact2] = Exact_sol_test(var_muu,Tipo)
syms( {'x', 'y','t' ,'muu',  'Pimp' , 'Pihalfmp'} )

symu_exact1 = -y*sin((x.^2 + y.^2)*sin(t^2+ 1)) -5;
symu_exact2 = x*sin((x.^2 + y.^2)*sin(t^2 + 1)) +10;
p_exact = sin((x+y)+t);

symf1 = diff(symu_exact1,t)-muu*laplacian(symu_exact1,[x,y]) + diff(p_exact,x);
symf2 = diff(symu_exact2,t)-muu*laplacian(symu_exact2,[x,y]) + diff(p_exact,y);

symL1 = -muu*laplacian(symu_exact1,[x,y]) ;
symL2 = -muu*laplacian(symu_exact2,[x,y]) ;

p_x_aux = matlabFunction(diff(p_exact,x),'vars',[x y t]);
p_x = @(puntos,t) p_x_aux(puntos(:,1),puntos(:,2),t);

p_y_aux = matlabFunction(diff(p_exact,y),'vars',[x y t]);
p_y = @(puntos,t) p_y_aux(puntos(:,1),puntos(:,2),t);

pii = Tipo('pi');
piihalf = Tipo('pi/2');


f1_aux = matlabFunction(symf1,'vars',[x y t muu Pimp Pihalfmp]);
f1 = @(puntos,t) f1_aux(puntos(:,1),puntos(:,2),t,var_muu,pii,piihalf);

f2_aux = matlabFunction(symf2,'vars',[x y t muu Pimp Pihalfmp]);
f2 = @(puntos,t) f2_aux(puntos(:,1),puntos(:,2),t,var_muu,pii,piihalf);


L1_aux = matlabFunction(symL1,'vars',[x y t muu Pimp Pihalfmp]);
L1 = @(puntos,t) L1_aux(puntos(:,1),puntos(:,2),t,var_muu,pii,piihalf);

L2_aux = matlabFunction(symL2,'vars',[x y t muu Pimp Pihalfmp]);
L2 = @(puntos,t) L2_aux(puntos(:,1),puntos(:,2),t,var_muu,pii,piihalf);

u_exact1_aux = matlabFunction(symu_exact1,'vars',[x y t muu Pimp Pihalfmp]);
u_exact1 = @(puntos,t) u_exact1_aux(puntos(:,1),puntos(:,2),t,var_muu,pii,piihalf);

u_exact2_aux = matlabFunction(symu_exact2,'vars',[x y t muu Pimp Pihalfmp]);
u_exact2 = @(puntos,t) u_exact2_aux(puntos(:,1),puntos(:,2),t,var_muu,pii,piihalf);

ut_exact1_aux = matlabFunction(diff(symu_exact1,t),'vars',[x y t muu Pimp Pihalfmp]);
ut_exact1 = @(puntos,t) ut_exact1_aux(puntos(:,1),puntos(:,2),t,var_muu,pii,piihalf);

ut_exact2_aux = matlabFunction(diff(symu_exact2,t),'vars',[x y t muu Pimp Pihalfmp]);
ut_exact2 = @(puntos,t) ut_exact2_aux(puntos(:,1),puntos(:,2),t,var_muu,pii,piihalf);


end
