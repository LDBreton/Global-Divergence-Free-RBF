%Funcion para obtener:
%f1,f2,dirchi1,dirchi2,p_exact_x,p_exact_y
%Usnado la solucion exacta de Wendland

function [f1,f2,u_exact1,u_exact2,B1,B2,p_x,p_y] = Exact_sol_slip(var_muu,Nx,Ny,Tipo)

syms( {'x', 'y','t','nX','nY','muu',  'Pimp' , 'Pihalfmp'} )

Sigma = @(f) (muu)*(jacobian([f(1) f(2)],[x,y]) ...
                    + transpose(jacobian([f(1) f(2)],[x,y])));
                
Sigman = @(f)  Sigma(f)*[nX ; nY];
Bn = @(f) f(1)*nX + f(2)*nY;
SlipB = @(f)  [-nY,nX]*Sigman(f);
 
symu_exact1 = -y*sin((x.^2 + y.^2)*2*Pimp)*sin(t*3*Pimp) + sin(t*20*Pimp);
symu_exact2 = x*sin((x.^2 + y.^2)*2*Pimp)*sin(t*3*Pimp) + cos(t*20*Pimp);
p_exact = sin((x+y)+t);

% symu_exact1 = -y*sin((x.^2 + y.^2)*2*Pimp)*sin(t*2*Pimp);
% symu_exact2 = x*sin((x.^2 + y.^2)*2*Pimp)*sin(t*2*Pimp);
% p_exact = sin((x+y)+t);

f1 = diff(symu_exact1,t)-muu*laplacian(symu_exact1,[x,y]) + diff(p_exact,x);
f2 = diff(symu_exact2,t)-muu*laplacian(symu_exact2,[x,y]) + diff(p_exact,y);

symB1 = Bn([symu_exact1,symu_exact2]);
symB2 = SlipB([symu_exact1,symu_exact2]);

pii = Tipo('pi');
piihalf = Tipo('pi/2');


f1_aux = matlabFunction(f1,'vars',[x y t muu Pimp Pihalfmp]);
f1 = @(puntos,t) f1_aux(puntos(:,1),puntos(:,2),t,var_muu,pii,piihalf);

f2_aux = matlabFunction(f2,'vars',[x y t muu Pimp Pihalfmp]);
f2 = @(puntos,t) f2_aux(puntos(:,1),puntos(:,2),t,var_muu,pii,piihalf);

B1_aux = matlabFunction(symB1,'vars',[x y t nX nY, muu Pimp Pihalfmp]);
B1 = @(puntos,t) B1_aux(puntos(:,1),puntos(:,2),t,...
                        Nx(puntos(:,1),puntos(:,2)),...
                        Ny(puntos(:,1),puntos(:,2)),...
                        var_muu,pii,piihalf);

B2_aux = matlabFunction(symB2,'vars',[x y t nX nY, muu Pimp Pihalfmp]);
B2 = @(puntos,t) B2_aux(puntos(:,1),puntos(:,2),t,...
                        Nx(puntos(:,1),puntos(:,2)),...
                        Ny(puntos(:,1),puntos(:,2)),...
                        var_muu,pii,piihalf);


u_exact1_aux = matlabFunction(symu_exact1,'vars',[x y t muu Pimp Pihalfmp]);
u_exact1 = @(puntos,t) u_exact1_aux(puntos(:,1),puntos(:,2),t,var_muu,pii,piihalf);

u_exact2_aux = matlabFunction(symu_exact2,'vars',[x y t muu Pimp Pihalfmp]);
u_exact2 = @(puntos,t) u_exact2_aux(puntos(:,1),puntos(:,2),t,var_muu,pii,piihalf);

p_x_aux = matlabFunction(diff(p_exact,x),'vars',[x y t]);
p_x = @(puntos,t) p_x_aux(puntos(:,1),puntos(:,2),t);

p_y_aux = matlabFunction(diff(p_exact,y),'vars',[x y t]);
p_y = @(puntos,t) p_y_aux(puntos(:,1),puntos(:,2),t);


end
