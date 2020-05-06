%for matlab 2016 

%Librerias
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('../Libs/AdvanpixMCT-4.4.5.12698') %% libreria de multipresicion
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

L1c = @(f) (-muu*laplacian(f(1),[x2,y2]) +diff(f(3),'x2'));
L2c = @(f) (-muu*laplacian(f(2),[x2,y2]) +diff(f(3),'y2'));

Ad1c = @(f) w1*diff(f(1),'x2') + w2*diff(f(1),'y2');
Ad2c = @(f) w1*diff(f(2),'x2') + w2*diff(f(2),'y2');

Ad1 = @(f) v1*diff(f(1),'x1') + v2*diff(f(1),'y1');
Ad2 = @(f) v1*diff(f(2),'x1') + v2*diff(f(2),'y1');

B1 = @(f) (f(1));
B2 = @(f) (f(2));
L1 = @(f) (-muu*laplacian(f(1),[x1,y1]) +diff(f(3),'x1'));
L2 = @(f) (-muu*laplacian(f(2),[x1,y1]) +diff(f(3),'y1'));

%%%contruyendo el Anzatz 
OperadoresY = {L1c,L2c,B1,B2,Ad1c,Ad2c}; %operadores para el anzatz
OperadoresX = {L1,L2,B1,B2,Ad1,Ad2}; %operadores que le vamos aplicar al anzatz
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


%Construyendo las matrices del sistema %SECCCION 3.2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[fbrGram,fbrAnzatz] = NavAnzat(fbrGram,fbrAnzatz);
U_zero = @(x,y) x.*0 + y.*0;
%Parte constante del sistema
Centros = {P_ni ,P_ni, P_f ,P_f };
M_phi1_const 	= BuildMGram(fbrAnzatz(1:4,1),{P_ni},Centros,Tipo,U_zero,U_zero);
M_phi2_const 	= BuildMGram(fbrAnzatz(1:4,2),{P_ni},Centros,Tipo,U_zero,U_zero);
M_Lphi1_const 	= BuildMGram(fbrGram(1:4,1),{P_ni},Centros,Tipo,U_zero,U_zero);
M_Lphi2_const 	= BuildMGram(fbrGram(1:4,2),{P_ni},Centros,Tipo,U_zero,U_zero);
M_Bphi1_const 	= BuildMGram(fbrGram(1:4,3),{P_f},Centros,Tipo,U_zero,U_zero);
M_Bphi2_const 	= BuildMGram(fbrGram(1:4,4),{P_f},Centros,Tipo,U_zero,U_zero);

%Parte no constante del sistema
M_phiAd1  =@(U1,U2) BuildMGram(fbrAnzatz(5:6,1),{P_ni},{P_ni,P_ni},Tipo,U1,U2);
M_phiAd2  =@(U1,U2) BuildMGram(fbrAnzatz(5:6,2),{P_ni},{P_ni,P_ni},Tipo,U1,U2);


M_LphiAd1 =@(U1,U2) BuildMGram(fbrGram(5:6,1),{P_ni},{P_ni,P_ni},Tipo,U1,U2) + ...
                    BuildMGram(fbrGram(5:6,5),{P_ni},{P_ni,P_ni},Tipo,U1,U2);
                
M_LphiAd2 =@(U1,U2) BuildMGram(fbrGram(5:6,2),{P_ni},{P_ni,P_ni},Tipo,U1,U2) + ...
                    BuildMGram(fbrGram(5:6,6),{P_ni},{P_ni,P_ni},Tipo,U1,U2);
                
M_BphiAd1 =@(U1,U2) BuildMGram(fbrGram(5:6,3),{P_f},{P_ni,P_ni},Tipo,U1,U2);
M_BphiAd2 =@(U1,U2) BuildMGram(fbrGram(5:6,4),{P_f},{P_ni,P_ni},Tipo,U1,U2);

Zero_aux1 = zeros(length(P_ni),2*length(P_f));
Zero_aux2 = zeros(length(P_f),2*length(P_f));


M_phi1 =@(U1,U2) M_phi1_const + [M_phiAd1(U1,U2),Zero_aux1];
M_phi2 =@(U1,U2) M_phi1_const + [M_phiAd2(U1,U2),Zero_aux1];
M_Lphi1 =@(U1,U2) M_Lphi1_const + [M_LphiAd1(U1,U2),Zero_aux1];
M_Lphi2 =@(U1,U2) M_Lphi2_const + [M_LphiAd2(U1,U2),Zero_aux1];
M_Bphi1 =@(U1,U2) M_Bphi1_const + [M_BphiAd1(U1,U2),Zero_aux2];
M_Bphi2 =@(U1,U2) M_Bphi2_const + [M_BphiAd2(U1,U2),Zero_aux2];


sistema_bdf1dt = @(dt,U1,U2) [M_phi1(U1,U2) + dt*M_Lphi1(U1,U2) ; ...
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
    [U1,U2] = interpaux(P_tot,u_exact1(P_tot,T),u_exact2(P_tot,T));
    %incrementando el tiempo
    b = RHS_bdf1(dt,T,yinicial1,yinicial2);
    M = sistema_bdf1dt(dt,U1,U2);
    %resolviendo el sistema
    Lambda = M\b;
    
    %obteniendo la solucion
    yinicial1= (M_phi1(U1,U2)*Lambda);
    yinicial2= (M_phi2(U1,U2)*Lambda);
    
    %calculando el error
    y1_exact = u_exact1(P_ni,T);
    y2_exact = u_exact2(P_ni,T);
    
    error1=max(abs(yinicial1-y1_exact));
    error2=max(abs(yinicial2-y2_exact));
    
    fprintf('El error del campo vector es %d \n',double(max(error1,error2)))
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%			   
