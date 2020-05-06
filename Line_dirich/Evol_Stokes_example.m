%% Librery
addpath('../Libs/AdvanpixMCT-4.5.2.12841') %% libreria de multipresicion
addpath('../Libs/DivFree_libreria/') %% libreria de las funciones
addpath('FreeFem_meshing/') %% libreria de la mallas

Loaddir = 'ConstRFBRStokes_wolfram/';
[fbrGram,fbrAnzatz] = LoadFbrFiles(Loaddir,4,8);

%% Tipo de precision en las variables
class_t = 'mp'; 
Tipo = @(foo) numeric_t(foo,class_t);
mp.Digits(100);


%% Obteniendo puntos de la malla
[ P_ni, P_f ] = Mesh_gen(1/7,1/7,1.0);
%%

%%
%transformamos los puntos a tipo Multi presicion
P_ni = Tipo(P_ni); % P_in = puntos interiores
P_f = Tipo(P_f);% P_f Puntos de Frontera


%% Construyendo las matrices del sistema
% Ver seccion 3.2 del articulo

%%
% Definiendo los centros las funciones de base radial
Centros = {P_ni ,P_ni, P_f ,P_f };

%%
% Se uliliza la funcion BuildMGram para construir las matrices
M_phi1 = BuildMGram(fbrAnzatz(:,1),{P_ni},Centros,Tipo);
M_phi2 = BuildMGram(fbrAnzatz(:,2),{P_ni},Centros,Tipo);
M_Lphi1 = BuildMGram(fbrGram(:,1),{P_ni},Centros,Tipo);
M_Lphi2 = BuildMGram(fbrGram(:,2),{P_ni},Centros,Tipo);
M_Bphi1 = BuildMGram(fbrGram(:,3),{P_f},Centros,Tipo);
M_Bphi2 = BuildMGram(fbrGram(:,4),{P_f},Centros,Tipo);
Interpolation_matrixPx = BuildMGram(fbrGram(:,5),{P_ni},Centros,Tipo);
Interpolation_matrixPy = BuildMGram(fbrGram(:,6),{P_ni},Centros,Tipo);

%% Programando las constantes del metodo BDF1 FBR
% Creando el sistema del lado izquiero del methodo bdf1
% en funcion del paso dt
sistema_bdf1dt = @(dt) [M_phi1 + dt*M_Lphi1 ; ...
                        M_phi2 + dt*M_Lphi2 ; ...
                        M_Bphi1; ...
                        M_Bphi2; ...
                        ];
%%
% Obteniendo la solucion exacta y las fuentes de la ecuacion
[f1,f2,u_exact1,u_exact2,p_x,p_y] = Exact_sol(var_muu,Tipo);

%%
% Creando el lado derecho del methodo bdf1
RHS_bdf1 = @(dt,T,Y1,Y2) vertcat( dt*f1(P_ni,T) + Y1, ...
								  dt*f2(P_ni,T) + Y2, ...
								  u_exact1(P_f,T),...
								  u_exact2(P_f,T));

%%
% Definiendo las constantes
steps = 5;
dt = mp('1/100');
Sistema_bdf1 = sistema_bdf1dt(dt);

%%
% Definiendo el tiempo inicial y las condiciones inicales
T = 0;
yinicial1= u_exact1(P_ni,T);
yinicial2= u_exact2(P_ni,T);

%% Aplicando el metodo BDF1
for i=1:steps
    disp(i);
    %incrementando el tiempo
    T = T + dt;
    b = RHS_bdf1(dt,T,yinicial1,yinicial2);
    
    %resolviendo el sistema
    Lambda = Sistema_bdf1\b;
    
    %obteniendo la solucion
    yinicial1= (M_phi1*Lambda);
    yinicial2= (M_phi2*Lambda);
    
    %calculando el error
    Errorpx= (Interpolation_matrixPx*Lambda)-p_x(P_ni,T);
    Errorpy= (Interpolation_matrixPy*Lambda)-p_y(P_ni,T);
    y1_exact = u_exact1(P_ni,T);
    y2_exact = u_exact2(P_ni,T);
    
    error1=max(abs(yinicial1-y1_exact));
    error2=max(abs(yinicial2-y2_exact));
    errorppx=max(abs(Errorpx));
    errorppy=max(abs(Errorpy));
    
    fprintf('El error del campo vector es %d \n',double(max(error1,error2)))
	fprintf('El de la precion es %d \n',double(max(errorppx,errorppy)))
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%			   
