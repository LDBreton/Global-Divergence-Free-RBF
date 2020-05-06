%for matlab 2016 

%Librerias
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('../Libs/AdvanpixMCT-4.5.2.12841') %% libreria de multipresicion
addpath('../Libs/DivFree_libreria/') %% libreria de las funciones
addpath('../Libs/distmesh/') %% libreria de la mallas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Definiendo el tipo de las variables para el calculo
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
class_t = 'mp'; 
Tipo = @(foo) numeric_t(foo,class_t);
mp.Digits(50);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%NOTA LA PARTE SYMBOLICA SE PUEDE CORRER UNA SOLA VEZ POR APARTE
%Parte Symbolica
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sym_vars = {'x1','y1','x2','y2','w1','w2','v1','v2','c','muu'};
var_c = Tipo('0.1');
var_muu = Tipo('1');
% centros = x2,y2
% nodos = x1,y1
syms(sym_vars); 
fbr = 1/sqrt(1 + c*((x1-x2)^2 + (y1-y2)^2));

disp('comenzando el calculo symbolico')
%%%operadores del Anzats%%%%%
%Operadores diferenciales del problema

L1c = @(f) (-muu*laplacian(f(1),[x2,y2]) +w1*diff(f(1),'x2') + w2*diff(f(1),'y2') +diff(f(3),'x2'));
L2c = @(f) (-muu*laplacian(f(2),[x2,y2]) +w1*diff(f(2),'x2') + w2*diff(f(2),'y2') +diff(f(3),'y2'));
B1 = @(f) (f(1));
B2 = @(f) (f(2));
L1 = @(f) (-muu*laplacian(f(1),[x1,y1])+ v1*diff(f(1),'x1') + v2*diff(f(1),'y1') +diff(f(3),'x1'));
L2 = @(f) (-muu*laplacian(f(2),[x1,y1])+ v1*diff(f(2),'x1') + v2*diff(f(2),'y1') +diff(f(3),'y1'));

%%%contruyendo el Anzatz 
OperadoresY = {L1c,L2c,B1,B2}; %operadores para el anzatz
OperadoresX = {L1,L2,B1,B2}; %operadores que le vamos aplicar al anzatz
[fbrGram,fbrAnzatz] = Fbr_Sym_Anzat(sym_vars,fbr,OperadoresY,OperadoresX,var_c,var_muu);
disp('fin del calculo symbolico')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Obteniendo puntos de la malla
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ P_ni, P_f ] = gridmakercirlceuni(10); 
%transformamos los puntos a tipo Multi presicion
P_ni = Tipo(P_ni); % P_in = puntos interiores
P_f = Tipo(P_f);% P_f Puntos de Frontera
P_tot = [P_ni;P_f];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[fbrGram,fbrAnzatz] = NavAnzat(fbrGram,fbrAnzatz);
U_zero = @(x,y) x.*0 + y.*0;
%Parte constante del sistema
Centros = {P_ni ,P_ni, P_f ,P_f };
M_phi1 	=@(U1,U2) BuildMGram(fbrAnzatz(:,1),{P_ni},Centros,Tipo,U1,U2);
M_phi2 	=@(U1,U2) BuildMGram(fbrAnzatz(:,2),{P_ni},Centros,Tipo,U1,U2);
M_Lphi1 =@(U1,U2) BuildMGram(fbrGram(:,1),{P_ni},Centros,Tipo,U1,U2);
M_Lphi2 =@(U1,U2) BuildMGram(fbrGram(:,2),{P_ni},Centros,Tipo,U1,U2);
M_Bphi1 =@(U1,U2) BuildMGram(fbrGram(:,3),{P_f},Centros,Tipo,U1,U2);
M_Bphi2 =@(U1,U2) BuildMGram(fbrGram(:,4),{P_f},Centros,Tipo,U1,U2);



sistema_bdf1 = @(dt,U1,U2) [M_phi1(U1,U2) + dt*M_Lphi1(U1,U2) ; ...
                              M_phi2(U1,U2) + dt*M_Lphi2(U1,U2) ; ...
                              M_Bphi1(U1,U2); ...
                              M_Bphi2(U1,U2); ...
                              ];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%Resolviendo y evolucionando el sistema
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[f1,f2,u_exact1,u_exact2,p_x,p_y] = Exact_sol(var_muu,Tipo);

RHS_bdf1 = @(dt,T,Y1,Y2) vertcat( dt*f1(P_ni,T) + Y1, ...
								  dt*f2(P_ni,T) + Y2, ...
								  u_exact1(P_f,T),...
								  u_exact2(P_f,T));

steps = 100;
dt = mp('1/100');

T = 0;
yinicial1= u_exact1(P_ni,T);
yinicial2= u_exact2(P_ni,T);

for i=1:steps
    disp(i);
    T = T + dt;
    b = RHS_bdf1(dt,T,yinicial1,yinicial2);
    
    for j=1:2
    %interpolando el campo de velocidades
    [U1,U2] = interpaux(P_tot,[yinicial1;u_exact1(P_f,T)],[yinicial2;u_exact2(P_f,T)]);
    %incrementando el tiempo
    M = sistema_bdf1(dt,U1,U2);
    %resolviendo el sistema
    Lambda = M\b;
    
    %obteniendo la solucion
    yinicial1= (M_phi1(U1,U2)*Lambda);
    yinicial2= (M_phi2(U1,U2)*Lambda);
    end
    %calculando el error
    y1_exact = u_exact1(P_ni,T);
    y2_exact = u_exact2(P_ni,T);
    
    error1=max(abs(yinicial1-y1_exact));
    error2=max(abs(yinicial2-y2_exact));
    
    fprintf('El error del campo vector es %d \n',double(max(error1,error2)))
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%			   
