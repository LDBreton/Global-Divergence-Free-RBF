% Librerias
 addpath('../Libs/AdvanpixMCT-4.5.2.12841') %% libreria de multipresicion
 addpath('../Libs/DivFree_libreria/');
 addpath('FreeFem_meshing/');
 SaveDire = 'FBRsaveFiles/';
[fbrGram,fbrAnzatz] = LoadFbrFiles(SaveDire,4,6);

class_t = 'mp'; 
Tipo = @(foo) numeric_t(foo,class_t);
mp.Digits(50);

var_muu = mp('1');
var_dt = mp('1/20');
var_c = mp('5');
Apoint = 1/20;
[P_ni,P_f] = Mesh_gen_N(50,0.5);
fprintf('%e %d \n',length(P_ni)/length(P_f),length(P_ni)+length(P_f));

% [mIdx,mD] = knnsearch([P_ni;P_f],[P_ni;P_f],'K',2);
% fprintf('Fill distance %e numero de puntos %d \n',max(mD(:,2)),length([P_ni;P_f]));

P_ni = Tipo(P_ni); % P_in = puntos interiores
P_f = Tipo(P_f);% P_f Puntos de Frontera


Centros = {P_ni ,P_ni, P_f ,P_f };


[f1,f2,u_exact1,u_exact2,p_x,p_y,L1,L2] = Exact_sol_test(var_muu,Tipo);

U_extact = @(t) [u_exact1(P_ni,t);u_exact2(P_ni,t)];
FU = @(t) [f1(P_ni,t);f2(P_ni,t)];
LU = @(t) [L1(P_ni,t);L2(P_ni,t)];
BU = @(t) [u_exact1(P_f,t);u_exact2(P_f,t)];

T = 0.0;
theta =mp('1');
Udiif = U_extact(T+var_dt) -U_extact(T);
URHS = -theta*var_dt*LU(T+var_dt)-(1-theta)*var_dt*LU(T) ...
       + var_dt*(theta*FU(T+var_dt) + (1-theta)*FU(T) );
    
fprintf('El Elocal es del campo velocidad %e \n',double(max(abs(Udiif-URHS))))
                              
                              
% Definiendo las constantes
Sistema_bdf1 = BuildMGram(fbrGram(1:4,1:4),Centros,Centros,Tipo,var_dt,var_c,var_muu);
M_phi1 = BuildMGram(fbrAnzatz(:,1),{P_ni},Centros,Tipo,var_dt,var_c,var_muu);
M_phi2 = BuildMGram(fbrAnzatz(:,2),{P_ni},Centros,Tipo,var_dt,var_c,var_muu);
Interpolation_matrixPx = BuildMGram(fbrGram(1:4,5),{P_ni},Centros,Tipo,var_dt,var_c,var_muu);
Interpolation_matrixPy = BuildMGram(fbrGram(1:4,6),{P_ni},Centros,Tipo,var_dt,var_c,var_muu);
% Definiendo el tiempo inicial y las condiciones inicales


T = 0;
yinicial1= u_exact1(P_ni,T);
yinicial2= u_exact2(P_ni,T);
steps = 1.0;
EvectY = double(zeros(steps,1));
EGradP = double(zeros(steps,1));
% Aplicando el metodo BDF1
for i=1:steps
%     disp(i);
    %incrementando el tiempo
    T = T + var_dt;
    b = [[yinicial1;yinicial2]+var_dt*FU(T) ; BU(T) ];
    
    %resolviendo el sistema
    Lambda = Sistema_bdf1\b;
    %obteniendo la solucion
    yinicial1= (M_phi1*Lambda);
    yinicial2= (M_phi2*Lambda);
    
    y1_exact = u_exact1(P_ni,T);
    y2_exact = u_exact2(P_ni,T);
    
    %calculando el error
    Errorpx= double((Interpolation_matrixPx*Lambda)-p_x(P_ni,T));
    Errorpy= double((Interpolation_matrixPy*Lambda)-p_y(P_ni,T));
 
    
    error1=max(abs(yinicial1-y1_exact));
    error2=max(abs(yinicial2-y2_exact));
    
    errorppx=max(abs(Errorpx));
    errorppy=max(abs(Errorpy));
    EvectY(i) = double(max(error1,error2));
    EGradP(i) = double(max(errorppx,errorppy));
    fprintf('El error del campo vector es %e \n',EvectY(i))
    fprintf('El de la precion es %e \n',EGradP(i))
end

pointsize = 30;
P_ni = double(P_ni); % P_in = puntos interiores
P_f = double(P_f);% P_f Puntos de Frontera             
             hold on
             plot(P_ni(:,1),P_ni(:,2),'.r')
             plot(P_f(:,1),P_f(:,2),'.b')
             scatter(P_ni(:,1), P_ni(:,2),pointsize,abs(Errorpx) + abs(Errorpy),'filled');
             %scatter(P_ni(:,1), P_ni(:,2),pointsize,abs(yinicial1-y1_exact) + abs(yinicial2-y2_exact),'filled');
             colorbar
             hold off
% pointsize = 10;
% scatter(P_ni(:,1), P_ni(:,2),pointsize,log(1+ abs(Errorpx)));
% 
% TRI = delaunay(double(P_ni(:,1)),double(P_ni(:,2)));
% trisurf(TRI,double(P_ni(:,1)),double(P_ni(:,2)),double(abs(Errorpx)))
% 
% hold on
% scatter3(P_ni(:,1), P_ni(:,2), abs(Errorpy));
% scatter3(P_ni(:,1), P_ni(:,2), abs(Errorpy));
% hold off