 addpath('../Libs/AdvanpixMCT-4.5.2.12841') %% libreria de multipresicion
 addpath('../Libs/DivFree_libreria/');
 addpath('FreeFem_meshing/');
 SaveDire = 'FBRsaveFilest/';
 mkdir(SaveDire)
% Tipo de precision en las variables
class_t = 'mp'; 
Tipo = @(foo) numeric_t(foo,class_t);
mp.Digits(50);

% Definiendo Variables simbolicas
% *Nota la parte sumbolica puede ser ejecutada una sola vez*
% centros = x2,y2 
% nodos = x1,y1
sym_vars = {'x1','y1','x2','y2','dtt','c','muu'};
assume(sym(sym_vars(1:end)),'real')
syms(sym_vars); 
%
% Las variables x1,y1,x2,y2 siempre tiene que estar al principio y
% presentes

% Definiendo la funcion de base radial
fbr = 1/sqrt(1 + c*((x1-x2)^2 + (y1-y2)^2));
%

% Definiendo los operadores de la ecuacion de Stokes 
% $$L_{1}(u,p) =  -\mu \Delta e_{1} + \frac{\partial e_3}{\partial x}$$
% 
% $$L_{2}(u,p) =  -\mu \Delta e_{2} + \frac{\partial e_3}{\partial y}$$
% 
% $$B_{1}(u,p) = e_{1}$$
% 
% $$B_{2}(u,p) = e_{2}$$
L1c = @(f) (-muu*laplacian(f(1),[x2,y2]) +diff(f(3),'x2'));
L2c = @(f) (-muu*laplacian(f(2),[x2,y2]) +diff(f(3),'y2'));
B1 = @(f) (f(1));
B2 = @(f) (f(2));
L1 = @(f) (f(1)-dtt*muu*laplacian(f(1),[x1,y1]) +dtt*diff(f(3),'x1'));
L2 = @(f) (f(2)-dtt*muu*laplacian(f(2),[x1,y1]) +dtt*diff(f(3),'y1'));

%
% Operadores para obtener la presion
% para obtener la presion
Px = @(f) diff(f(3),'x1');
Py = @(f) diff(f(3),'y1');
%

%
% Juntando los operadores en 2 celdas
OperadoresY = {L1c,L2c,B1,B2}; %operadores para el anzatz
OperadoresX = {L1,L2,B1,B2,Px,Py}; %operadores que le vamos aplicar al anzatz

% Realizando el calculo simbolico
disp('comenzando el calculo symbolico')
[fbrGram,fbrAnzatz] = SymBuildMGram(SaveDire,sym_vars,fbr,OperadoresY,OperadoresX);
disp('fin del calculo symbolico')