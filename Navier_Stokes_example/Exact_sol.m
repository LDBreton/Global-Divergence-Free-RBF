%Funcion para obtener:
%f1,f2,dirchi1,dirchi2,p_exact_x,p_exact_y
%Usnado la solucion exacta de Wendland

function [f1,f2,u_exact1,u_exact2,p_x,p_y] = Exact_sol(var_muu,Tipo)
syms( {'x', 'y' , 't' , 'mus',  'Pimp' , 'Pihalfmp'} )


Symu_exact1 = -Pimp*y*cos(Pihalfmp*(x.^2 + y.^2))*cos(Pimp*t);
Symu_exact2 = Pimp*x*cos(Pihalfmp*(x.^2 + y.^2))*cos(Pimp*t);
p_exact = sin((x-y)+t);
%p_exact = x+y;
f1 = diff(Symu_exact1,t)-mus*laplacian(Symu_exact1,[x,y]) + Symu_exact1*diff(Symu_exact1,'x') +Symu_exact2*diff(Symu_exact1,'y') + diff(p_exact,x);
f2 = diff(Symu_exact2,t)-mus*laplacian(Symu_exact2,[x,y]) +  Symu_exact1*diff(Symu_exact2,'x') +Symu_exact2*diff(Symu_exact2,'y') + diff(p_exact,y);

p_x_aux = matlabFunction(diff(p_exact,x),'vars',[x y t]);
p_x = @(puntos,t) p_x_aux(puntos(:,1),puntos(:,2),t);

p_y_aux = matlabFunction(diff(p_exact,y),'vars',[x y t]);
p_y = @(puntos,t) p_y_aux(puntos(:,1),puntos(:,2),t);

pii = Tipo('pi');
piihalf = Tipo('pi/2');


f1_aux = matlabFunction(f1,'vars',[x y t mus Pimp Pihalfmp]);
f1 = @(puntos,t) f1_aux(puntos(:,1),puntos(:,2),t,var_muu,pii,piihalf);

f2_aux = matlabFunction(f2,'vars',[x y t mus Pimp Pihalfmp]);
f2 = @(puntos,t) f2_aux(puntos(:,1),puntos(:,2),t,var_muu,pii,piihalf);

u_exact1_aux = matlabFunction(Symu_exact1,'vars',[x y t mus Pimp Pihalfmp]);
u_exact1 = @(puntos,t) u_exact1_aux(puntos(:,1),puntos(:,2),t,var_muu,pii,piihalf);

u_exact2_aux = matlabFunction(Symu_exact2,'vars',[x y t mus Pimp Pihalfmp]);
u_exact2 = @(puntos,t) u_exact2_aux(puntos(:,1),puntos(:,2),t,var_muu,pii,piihalf);


end
