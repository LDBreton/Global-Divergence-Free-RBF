%% Librerias
 addpath('../Libs/AdvanpixMCT-4.5.2.12841') %% libreria de multipresicion
 addpath('../Libs/DivFree_libreria/');
 addpath('FreeFem_meshing/');
 SaveDire = 'FBRsaveFiles/';
[fbrGram,fbrAnzatz] = LoadFbrFiles(SaveDire,4,6);

class_t = 'mp'; 
Tipo = @(foo) numeric_t(foo,class_t);
mp.Digits(200);
%%

% Obteniendo la solucion exacta y las fuentes de la ecuacion
Avar_muu = mp('1.0');
Avar_dt = mp('1/10');
Apoint = 1/15;

for lk=1:length(Avar_muu)
    for k=1:length(Apoint)
        for j=1:length(Avar_dt)
            [P_ni,P_f] = Mesh_gen(Apoint(k)*0.4,Apoint(k),0.5);
             var_muu = Avar_muu(lk);
             var_dt = Avar_dt(j);
             var_c = mp('0.2');

            P_ni = Tipo(P_ni); % P_in = puntos interiores
            P_f = Tipo(P_f);% P_f Puntos de Frontera
            
            Centros = {P_ni ,P_ni, P_f ,P_f };
            
            [f1,f2,u_exact1,u_exact2,p_x,p_y] = Exact_sol(var_muu,Tipo);
            
            %%
            % Creando el lado derecho del methodo bdf1
            
            RHS_bdf1 = @(dt,T,Y1,Y2) vertcat( dt*f1(P_ni,T) + Y1, ...
                dt*f2(P_ni,T) + Y2, ...
                u_exact1(P_f,T),...
                u_exact2(P_f,T));
            
            
            %%
            % Definiendo las constantes
            Sistema_bdf1 = BuildMGram(fbrGram(1:4,1:4),Centros,Centros,Tipo,var_dt,var_c,var_muu);
            M_phi1 = BuildMGram(fbrAnzatz(:,1),{P_ni},Centros,Tipo,var_dt,var_c,var_muu);
            M_phi2 = BuildMGram(fbrAnzatz(:,2),{P_ni},Centros,Tipo,var_dt,var_c,var_muu);
            Interpolation_matrixPx = BuildMGram(fbrGram(1:4,5),{P_ni},Centros,Tipo,var_dt,var_c,var_muu);
            Interpolation_matrixPy = BuildMGram(fbrGram(1:4,6),{P_ni},Centros,Tipo,var_dt,var_c,var_muu);
            %%
            
            % Definiendo el tiempo inicial y las condiciones inicales
            T = 0;
            yinicial1= u_exact1(P_ni,T);
            yinicial2= u_exact2(P_ni,T);
            Lambda = zeros(2*length(P_ni) + 2*length(P_f),1);
            
            steps = 1.0;
            EvectY = double(zeros(steps,1));
            EGradP = double(zeros(steps,1));

            for i=1:steps
                disp(i);
                %incrementando el tiempo
                T = T + var_dt;
                b = RHS_bdf1(var_dt,T,yinicial1,yinicial2);
                Lambda = Sistema_bdf1\b;
                
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
                
                EvectY(i) = double(max(error1,error2));
                EGradP(i) = double(max(errorppx,errorppy));
                fprintf('El error del campo vector es %e \n',EvectY(i))
                fprintf('El de la precion es %e \n',EGradP(i))
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
    end
end