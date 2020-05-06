% Librerias
addpath('../Libs/AdvanpixMCT-4.5.2.12841') %% libreria de multipresicion
addpath('../Libs/DivFree_libreria/');
addpath('FreeFem_meshing/');
addpath('locallibs/')
Loaddir = 'FBRLines/';
[fbrGram,fbrAnzatz] = LoadFbrFilesV2(Loaddir,4,6);
Savetable = 'NewTables/';
mkdir(Savetable);
class_t = 'mp';
Tipo = @(foo) numeric_t(foo,class_t);
mp.Digits(100);


muu = Tipo('1');
ParamC = Tipo('0.1');

%%%%%%%%%%%%%%%%%%%%%%%%%%PARAMETROS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ P_ni, P_f ] = Mesh_gen(1/10,1/10,0.5);
 %crea una malla de NxN en el [0,1]

P_ni = Tipo(P_ni); % P_in = puntos interiores
P_f = Tipo((P_f));% P_f Puntos de Frontera
P_tot = [P_ni;P_f];
Centros = {P_ni ,P_ni, P_f ,P_f};


Interpolation_matrix1 = BuildMGram(fbrAnzatz(:,1),{P_ni},Centros,Tipo,ParamC,muu);
BInterpolation_matrix1 = BuildMGram(fbrAnzatz(:,1),{P_f},Centros,Tipo,ParamC,muu);
Interpolation_matrix2 = BuildMGram(fbrAnzatz(:,2),{P_ni},Centros,Tipo,ParamC,muu);
BInterpolation_matrix2 = BuildMGram(fbrAnzatz(:,2),{P_f},Centros,Tipo,ParamC,muu);


Interpolation_matrixPx = BuildMGram(fbrGram(:,5),{P_ni},Centros,Tipo,ParamC,muu);
Interpolation_matrixPy = BuildMGram(fbrGram(:,6),{P_ni},Centros,Tipo,ParamC,muu);

L1_matrix = BuildMGram(fbrGram(:,1),{P_ni},Centros,Tipo,ParamC,muu);
L2_matrix = BuildMGram(fbrGram(:,2),{P_ni},Centros,Tipo,ParamC,muu);


times_vect = [10,100,500,1000,5000,10000,20000]';

Error_infty_u = zeros(length(times_vect),1);
Error_grad_p = zeros(length(times_vect),1);
for ite = 1:length(times_vect)

dt = Tipo(strcat('1/',num2str(times_vect(ite))));


sistema = [Interpolation_matrix1 + dt*L1_matrix ; ...
           Interpolation_matrix2 + dt*L2_matrix ; ...
           BInterpolation_matrix1; ...
           BInterpolation_matrix2; ...
           ];   
sistema1 = [Interpolation_matrix1 + mp('2/3')*dt*L1_matrix ; ...
           Interpolation_matrix2 + mp('2/3')*dt*L2_matrix ; ...
           BInterpolation_matrix1; ...
           BInterpolation_matrix2; ...
           ];   
              
       
rhs1_matrix = Interpolation_matrix1 ;
rhs2_matrix = Interpolation_matrix2 ;
       
       
[f1,f2,dirchi1,dirchi2,p_x,p_y] = circle_exact_sol(muu);
yinicial = [dirchi1(P_ni(:,1),P_ni(:,2),0); ...
            dirchi2(P_ni(:,1),P_ni(:,2),0); ...
            dirchi1(P_f(:,1),P_f(:,2),0); ...
            dirchi2(P_f(:,1),P_f(:,2),0) ] ;

A = [Interpolation_matrix1;...
     Interpolation_matrix2;...
     BInterpolation_matrix1;...
     BInterpolation_matrix2;];
%double(cond(A))
lamda1 = A\yinicial;


%[U,S,V]=svd(sistema,0);
[U1,S1,V1]=svd(sistema1,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Fin de la creacion la matriz de Gram')



%%%%%%%%%%%%%%%%%%%%%%%%RESOLVIENDO EL SISTEMA%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Resolviendo el sistema')
steps = floor(1/double(dt));
times = Tipo(linspace(mp(0),mp(1),steps+1));
error = zeros(steps,1);
errorp = zeros(steps,1);

yinicial1 = dirchi1(P_ni(:,1),P_ni(:,2),0);
yinicial2 = dirchi2(P_ni(:,1),P_ni(:,2),0);
for i=1:1
    
    
    tiemp = Tipo(zeros(length(P_ni(:,1)),1))+times(i)+dt;
    
    y1_exact = dirchi1(P_ni(:,1),P_ni(:,2),tiemp);
    y2_exact = dirchi2(P_ni(:,1),P_ni(:,2),tiemp);
    u1 = f1(P_ni(:,1),P_ni(:,2),tiemp);
    u2 = f2(P_ni(:,1),P_ni(:,2),tiemp);
    
    tiemp = zeros(length(P_f(:,1)),1)+times(i)+dt;
    b1 = dirchi1(P_f(:,1),P_f(:,2),tiemp);
    b2 = dirchi2(P_f(:,1),P_f(:,2),tiemp);
    
    
    disp(i);
   b =  vertcat( dt*u1 + rhs1_matrix*lamda1, ...
                 dt*u2 + rhs2_matrix*lamda1, ...
                 b1,...
                 b2);    
 %lamda2= V*((U'*b)./diag(S));
  lamda2 = sistema\b;

yinicial1= (Interpolation_matrix1*lamda2);
yinicial2= (Interpolation_matrix2*lamda2);

Errorpx= (Interpolation_matrixPx*lamda2)-p_x(P_ni(:,1),P_ni(:,2),times(i)+dt);
Errorpy= (Interpolation_matrixPy*lamda2)-p_y(P_ni(:,1),P_ni(:,2),times(i)+dt);

error1=max(abs(yinicial1-y1_exact));
error2=max(abs(yinicial2-y2_exact));
errorppx=max(abs(Errorpx));
errorppy=max(abs(Errorpy));

error(i) = max(error1,error2);
errorp(i) = max(errorppx,errorppy);
fprintf('erroy %2.3e errop %2.3e \n',max(error1,error2),max(errorppx,errorppy));
end


for i=2:steps
    
    
    tiemp = Tipo(zeros(length(P_ni(:,1)),1))+times(i)+dt;
    
    y1_exact = dirchi1(P_ni(:,1),P_ni(:,2),tiemp);
    y2_exact = dirchi2(P_ni(:,1),P_ni(:,2),tiemp);
    u1 = f1(P_ni(:,1),P_ni(:,2),tiemp);
    u2 = f2(P_ni(:,1),P_ni(:,2),tiemp);
    
    tiemp = zeros(length(P_f(:,1)),1)+times(i)+dt;
    b1 = dirchi1(P_f(:,1),P_f(:,2),tiemp);
    b2 = dirchi2(P_f(:,1),P_f(:,2),tiemp);
    
    
   b= vertcat(mp('2/3')*dt*u1 +mp('4/3')*rhs1_matrix*lamda2 -mp('1/3')*rhs1_matrix*lamda1,...
              mp('2/3')*dt*u2 +mp('4/3')*rhs2_matrix*lamda2 -mp('1/3')*rhs2_matrix*lamda1,...
              b1,...
              b2);     
          
 lamda1 = lamda2;     
 lamda2= V1*((U1'*b)./diag(S1));
 %lamda2 = sistema1\b;

yinicial1= (Interpolation_matrix1*lamda2);
yinicial2= (Interpolation_matrix2*lamda2);

Errorpx= (Interpolation_matrixPx*lamda2)-p_x(P_ni(:,1),P_ni(:,2),times(i)+dt);
Errorpy= (Interpolation_matrixPy*lamda2)-p_y(P_ni(:,1),P_ni(:,2),times(i)+dt);

error1=max(abs(yinicial1-y1_exact));
error2=max(abs(yinicial2-y2_exact));
errorppx=max(abs(Errorpx));
errorppy=max(abs(Errorpy));

error(i) = max(error1,error2);
errorp(i) = max(errorppx,errorppy);
fprintf('erroy %2.3e errop %2.3e \n',max(error1,error2),max(errorppx,errorppy));
end


fprintf('Max erroy %2.3e Max errop %2.3e \n',...
    max(error),max(errorp));

Error_infty_u(ite) = double(max(error));
Error_grad_p(ite) = double(max(errorp));

erroresmat = [1./times_vect , Error_infty_u , Error_grad_p];
dlmwrite('errores_bdf2_mu_e0.txt',erroresmat,'delimiter','\t','precision',10);

end
