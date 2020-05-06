%Funcion para obtener:
%f1,f2,dirchi1,dirchi2,p_exact_x,p_exact_y
%Usnado la solucion exacta de Wendland

function [f1M,f2M,u_exact1M,u_exact2M,p_x,p_y] = Exact_sol2(var_muu,Tipo)
syms( {'x', 'y','t' ,'muu',  'Pimp' , 'Pihalfmp'} )

symu_exact1 = -y*sin((x.^2 + y.^2)*2*Pimp)*sin(t*2*Pimp);
symu_exact2 = x*sin((x.^2 + y.^2)*2*Pimp)*sin(t*2*Pimp);
p_exact = sin((x+y)+t);

f1 = diff(symu_exact1,t)-muu*laplacian(symu_exact1,[x,y]) + diff(p_exact,x);
f2 = diff(symu_exact2,t)-muu*laplacian(symu_exact2,[x,y]) + diff(p_exact,y);

p_x_aux = matlabFunction(diff(p_exact,x),'vars',[x y t]);
p_x = @(puntos,t) p_x_aux(puntos(:,1),puntos(:,2),t);

p_y_aux = matlabFunction(diff(p_exact,y),'vars',[x y t]);
p_y = @(puntos,t) p_y_aux(puntos(:,1),puntos(:,2),t);

pii = Tipo('pi');
piihalf = Tipo('pi/2');


f1_aux = matlabFunction(f1,'vars',[x y t muu Pimp Pihalfmp]);
f1M = @(puntos,t) f1_aux(puntos(:,1),puntos(:,2),t,var_muu,pii,piihalf);

f2_aux = matlabFunction(f2,'vars',[x y t muu Pimp Pihalfmp]);
f2M = @(puntos,t) f2_aux(puntos(:,1),puntos(:,2),t,var_muu,pii,piihalf);

u_exact1_aux = matlabFunction(symu_exact1,'vars',[x y t muu Pimp Pihalfmp]);
u_exact1M = @(puntos,t) u_exact1_aux(puntos(:,1),puntos(:,2),t,var_muu,pii,piihalf);

u_exact2_aux = matlabFunction(symu_exact2,'vars',[x y t muu Pimp Pihalfmp]);
u_exact2M = @(puntos,t) u_exact2_aux(puntos(:,1),puntos(:,2),t,var_muu,pii,piihalf);


end
