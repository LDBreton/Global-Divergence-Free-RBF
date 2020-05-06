%% Librerias
addpath('../Libs/AdvanpixMCT-4.5.2.12841') %% libreria de multipresicion
addpath('../Libs/DivFree_libreria/');
addpath('FreeFem_meshing/');
addpath('locallibs/')
Loaddir = 'FBRsaveFiles/';
[fbrGram,fbrAnzatz] = LoadFbrFiles(Loaddir,4,8);
Savetable = 'NewTables/';
mkdir(Savetable);
class_t = 'mp';
Tipo = @(foo) numeric_t(foo,class_t);
mp.Digits(100);
%%

% Obteniendo la solucion exacta y las fuentes de la ecuacion
Avar_muu = [mp('1.0'),mp('1e-06')];
%,mp('1e-06')]
Avar_dt = [mp('0.01'),mp('0.001'),mp('0.0001')];
Apoint = [1/5,1/10,1/15,1/17];
for lk=1:length(Avar_muu)
    fileID = fopen([Savetable,'Errores_muu_',num2str(double(Avar_muu(lk))),'_dirichbdf2_test.txt'],'w');
    for j=1:length(Avar_dt)
        for k=1:length(Apoint)
            [P_ni,P_f] = Mesh_gen(Apoint(k),Apoint(k),1.0);
            var_muu = Avar_muu(lk);
            var_dt = Avar_dt(j);
            var_c = mp('0.5');
            fprintf('%d %d ',length(P_ni),length(P_f));
            
            fprintf(fileID,'%d %5.5e ',double(length([P_ni;P_f])),double(var_dt));
            
            P_ni = Tipo(P_ni); % P_in = puntos interiores
            P_f = Tipo(P_f);% P_f Puntos de Frontera
            
            Centros = {P_ni ,P_ni, P_f ,P_f };
            
            [f1,f2,dirchi1,dirchi2,p_x,p_y] = Exact_sol2(var_muu,Tipo);
           
            
            
            %%
            % Creando el lado derecho del methodo bdf1
            
            % Definiendo las constantes
            Sistema_bdf2 = BuildMGram(fbrGram(1:4,1:4),Centros,Centros,Tipo,mp('2/3')*var_dt,var_c,var_muu);
            M_phi1_bdf2 = BuildMGram(fbrAnzatz(:,1),{P_ni},Centros,Tipo,mp('2/3')*var_dt,var_c,var_muu);
            M_phi2_bdf2 = BuildMGram(fbrAnzatz(:,2),{P_ni},Centros,Tipo,mp('2/3')*var_dt,var_c,var_muu);
            Interpolation_matrixPx_bdf2 = BuildMGram(fbrGram(1:4,5),{P_ni},Centros,Tipo,mp('2/3')*var_dt,var_c,var_muu);
            Interpolation_matrixPy_bdf2 = BuildMGram(fbrGram(1:4,6),{P_ni},Centros,Tipo,mp('2/3')*var_dt,var_c,var_muu);
            %%
            
            % Definiendo el tiempo inicial y las condiciones inicales
            tiemp = var_dt;
            yinicial_n_0_1= dirchi1(P_ni,tiemp);
            yinicial_n_0_2= dirchi2(P_ni,tiemp);
            yinicial_n_1_1= dirchi1(P_ni,tiemp+ var_dt);
            yinicial_n_1_2= dirchi2(P_ni,tiemp+ var_dt);
            tiemp = tiemp + var_dt;
            
            steps = floor(1/var_dt)-1;
            EvectY = double(zeros(steps,1));
            EGradP = double(zeros(steps,1));
            % Aplicando el metodo BDF1
            if(steps >= 99)
                [U,S,V]=svd(Sistema_bdf2);
            end
            for i=2:steps
                disp(i);
                tiemp = tiemp + var_dt;
                
                FF1 = f1(P_ni,tiemp);
                FF2 = f2(P_ni,tiemp);
                b1 = dirchi1(P_f,tiemp);
                b2 = dirchi2(P_f,tiemp);
                b= vertcat(mp('2./3.')*var_dt*FF1 +mp('4/3')*yinicial_n_1_1 -mp('1/3')*yinicial_n_0_1,...
                    mp('2./3.')*var_dt*FF2 +mp('4/3')*yinicial_n_1_2 -mp('1/3')*yinicial_n_0_2,...
                    b1,...
                    b2);
                
                if(steps >= 99)
                    Lambda= V*((U'*b)./diag(S));
                else
                    Lambda = Sistema_bdf2\b;
                end
                yinicial_n_0_1 = yinicial_n_1_1;
                yinicial_n_0_2 = yinicial_n_1_2;
                %incrementando el tiempo
                %resolviendo el sistema
                
                %obteniendo la solucion
                yinicial_n_1_1= (M_phi1_bdf2*Lambda);
                yinicial_n_1_2= (M_phi2_bdf2*Lambda);
                
                %calculando el error
                Errorpx= (Interpolation_matrixPx_bdf2*Lambda)-p_x(P_ni,tiemp);
                Errorpy= (Interpolation_matrixPy_bdf2*Lambda)-p_y(P_ni,tiemp);
                y1_exact = dirchi1(P_ni,tiemp);
                y2_exact = dirchi2(P_ni,tiemp);
                
                error1=max(abs(yinicial_n_1_1-y1_exact));
                error2=max(abs(yinicial_n_1_2-y2_exact));
                
                errorppx=max(abs(Errorpx));
                errorppy=max(abs(Errorpy));
                
                EvectY(i) = EvectY(i) + double(max(error1,error2));
                EGradP(i) = double(max(errorppx,errorppy));
                fprintf('El error del campo vector es %e \n',EvectY(i))
                fprintf('El de la precion es %e \n',EGradP(i))
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            fprintf(fileID,'%e %e %e \n',max(EvectY),max(EGradP),double(cond(Sistema_bdf2)));
        end
    end
    fclose(fileID);
end