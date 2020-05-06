addpath('AdvanpixMCT-4.4.5.12698')
addpath('DivFree_libreria/') %% libreria de las funciones
addpath('distmesh/')
mp.Digits(50);

%NOTA LA PARTE SYMBOLICA SE PUEDE CORRER UNA SOLA VEZ POR APARTE
%%%%%%%%%%%%%%%%%%%%%%%PARTE SYMBOLICA%%%%%%%%%%%%%%%%%%%%%
Tipo = @(x) mp(x); %%si quieremos cambiar a multipre
ParamC = Tipo('0.1'); %%% PARAMETRO DE FORMA DE LA FUNCION SI TIENE
N = 10; %calidad del mallado
%Variables symbolicas
sym_vars = {'x','y','c','mus'};
syms(sym_vars);

%fbr = exp(-(c^2)*(x^2 + y^2));
%fbr = -((x^2 + y^2))^6 * log((x^2 + y^2))/2;
fbr = 1/sqrt(1 + c*(x^2 + y^2));

dt = Tipo('1/100');
muu = Tipo('1');

disp('comenzando el calculo symbolico')
%Operadores diferenciales del problema
L1y = @(f) (-mus*laplacian(f(1),[x,y]) -diff(f(3),'x'));
L2y = @(f) (-mus*laplacian(f(2),[x,y]) -diff(f(3),'y'));
B1 = @(f) (f(1));
B2 = @(f) (f(2));
L1 = @(f) (-mus*laplacian(f(1),[x,y]) +diff(f(3),'x'));
L2 = @(f) (-mus*laplacian(f(2),[x,y]) +diff(f(3),'y'));
Px = @(f) diff(f(3),'x');
Py = @(f) diff(f(3),'y');


%%%contruyendo el Anzatz VER BuildAnzatzmp.M
OperadoresY = {L1y,L2y,B1,B2};
OperadoresX = {L1,L2,B1,B2,Px,Py};
[fbrGram,fbrAnzatz] = BuildAnzatzmp(sym_vars,fbr,OperadoresY,OperadoresX,ParamC,muu);
disp('fin del calculo symbolico')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%PARAMETROS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ P_ni, P_f ] = gridmakercirlceuni(N); %crea una malla de NxN en el [0,1]

P_ni = Tipo(P_ni); % P_in = puntos interiores
P_f = Tipo((P_f));% P_f Puntos de Frontera
P_tot = [P_ni;P_f];
Centros = {P_ni ,P_ni, P_f ,P_f};



disp('Creando las matrices del sistema')
%%%creando las matrices del articulo%%%%%%%
Interpolation_matrix1 = BuildMGram(fbrAnzatz(:,1),{P_ni},Centros,Tipo);
BInterpolation_matrix1 = BuildMGram(fbrAnzatz(:,1),{P_f},Centros,Tipo);
Interpolation_matrix2 = BuildMGram(fbrAnzatz(:,2),{P_ni},Centros,Tipo);
BInterpolation_matrix2 = BuildMGram(fbrAnzatz(:,2),{P_f},Centros,Tipo);

Interpolation_matrixPx = BuildMGram(fbrGram(:,5),{P_ni},Centros,Tipo);
Interpolation_matrixPy = BuildMGram(fbrGram(:,6),{P_ni},Centros,Tipo);

L1_matrix = BuildMGram(fbrGram(:,1),{P_ni},Centros,Tipo);
L2_matrix = BuildMGram(fbrGram(:,2),{P_ni},Centros,Tipo);

sistema_bdf1 = [Interpolation_matrix1 + dt*L1_matrix ; ...
    Interpolation_matrix2 + dt*L2_matrix ; ...
    BInterpolation_matrix1; ...
    BInterpolation_matrix2; ...
    ];

sistema_bdf2 = [Interpolation_matrix1 + mp('2/3')*dt*L1_matrix ; ...
    Interpolation_matrix2  + mp('2/3')*dt*L2_matrix ; ...
    BInterpolation_matrix1; ...
    BInterpolation_matrix2; ...
    ];


rhs1_matrix = Interpolation_matrix1 ;
rhs2_matrix = Interpolation_matrix2 ;
%%%%%%%%%%%%%%%%%%%%%%%%%
disp('fin de la creacion de las matrices del sistema')


%Obteniendo fuente y condiciones de frontera de la solucion
[f1,f2,dirchi1,dirchi2,p_x,p_y] = circle_exact_sol(muu);

%obteniendo el primer coefficiente lamda(0)
yinicial = [dirchi1(P_ni(:,1),P_ni(:,2),0); ...
    dirchi2(P_ni(:,1),P_ni(:,2),0); ...
    dirchi1(P_f(:,1),P_f(:,2),0); ...
    dirchi2(P_f(:,1),P_f(:,2),0) ] ;
A = [Interpolation_matrix1;...
    Interpolation_matrix2;...
    BInterpolation_matrix1;...
    BInterpolation_matrix2;];
lamda1 = A\yinicial;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%RESOLVIENDO EL SISTEMA%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Resolviendo el sistema')
steps = floor(1/double(dt));
times = Tipo(linspace(mp(0),mp(1),steps+1));
error = zeros(steps,1);
errorp = zeros(steps,1);

yinicial1 = dirchi1(P_ni(:,1),P_ni(:,2),0);
yinicial2 = dirchi2(P_ni(:,1),P_ni(:,2),0);
for i=1:1
    
    
    tiemp = times(i)+dt;
    u1 = f1(P_ni(:,1),P_ni(:,2),tiemp);
    u2 = f2(P_ni(:,1),P_ni(:,2),tiemp);
    b1 = dirchi1(P_f(:,1),P_f(:,2),tiemp);
    b2 = dirchi2(P_f(:,1),P_f(:,2),tiemp);
    
    
    disp(i);
    b =  vertcat( dt*u1 + rhs1_matrix*lamda1, ...
        dt*u2 + rhs2_matrix*lamda1, ...
        b1,...
        b2);
    
    %resolviendo el sistema
    lamda2 = sistema_bdf1\b;
    
    %obteniendo la solucion
    yinicial1= (Interpolation_matrix1*lamda2);
    yinicial2= (Interpolation_matrix2*lamda2);
    
    Errorpx= (Interpolation_matrixPx*lamda2)-p_x(P_ni(:,1),P_ni(:,2),times(i)+dt);
    Errorpy= (Interpolation_matrixPy*lamda2)-p_y(P_ni(:,1),P_ni(:,2),times(i)+dt);
    
    y1_exact = dirchi1(P_ni(:,1),P_ni(:,2),tiemp);
    y2_exact = dirchi2(P_ni(:,1),P_ni(:,2),tiemp);
    
    error1=max(abs(yinicial1-y1_exact));
    error2=max(abs(yinicial2-y2_exact));
    errorppx=max(abs(Errorpx));
    errorppy=max(abs(Errorpy));
    
    max(error1,error2)
    max(errorppx,errorppy)
    error(i) = max(error1,error2);
    errorp(i) = max(errorppx,errorppy);
end


for i=2:steps
    
    
    tiemp = times(i)+dt;
    
    u1 = f1(P_ni(:,1),P_ni(:,2),tiemp);
    u2 = f2(P_ni(:,1),P_ni(:,2),tiemp);
    b1 = dirchi1(P_f(:,1),P_f(:,2),tiemp);
    b2 = dirchi2(P_f(:,1),P_f(:,2),tiemp);
    
    
    b= vertcat(mp('2/3')*dt*u1 +mp('4/3')*dirchi1(P_ni(:,1),P_ni(:,2),times(i)) -mp('1/3')*dirchi1(P_ni(:,1),P_ni(:,2),times(i-1)),...
        mp('2/3')*dt*u2 +mp('4/3')*dirchi2(P_ni(:,1),P_ni(:,2),times(i)) -mp('1/3')*dirchi2(P_ni(:,1),P_ni(:,2),times(i-1)),...
        b1,...
        b2);
    
    lamda1 = lamda2;
    
    %resolviendo el sistema
    lamda2 = sistema_bdf2\b;
    
    %obteniendo la solucion
    yinicial1= (Interpolation_matrix1*lamda2);
    yinicial2= (Interpolation_matrix2*lamda2);
    
    Errorpx= (Interpolation_matrixPx*lamda2)-p_x(P_ni(:,1),P_ni(:,2),times(i)+dt);
    Errorpy= (Interpolation_matrixPy*lamda2)-p_y(P_ni(:,1),P_ni(:,2),times(i)+dt);
    
    
    y1_exact = dirchi1(P_ni(:,1),P_ni(:,2),tiemp);
    y2_exact = dirchi2(P_ni(:,1),P_ni(:,2),tiemp);
    
    error1=max(abs(yinicial1-y1_exact));
    error2=max(abs(yinicial2-y2_exact));
    errorppx=max(abs(Errorpx));
    errorppy=max(abs(Errorpy));
    
    error(i) = max(error1,error2);
    errorp(i) = max(errorppx,errorppy);
    
    
    fprintf('El error del campo vector es %d \n',double(max(error1,error2)))
	fprintf('El de la precion es %d \n',double(max(errorppx,errorppy)))
end
