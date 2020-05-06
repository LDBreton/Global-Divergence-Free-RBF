 addpath('../Libs/AdvanpixMCT-4.5.2.12841') %% libreria de multipresicion
 addpath('../Libs/DivFree_libreria/');
 addpath('FreeFem_meshing/');
 SaveDire = 'FBRsave_wolfram/';
 mkdir(SaveDire)
 %Tipo de precision en las variables
 class_t = 'mp'; 
 Tipo = @(foo) numeric_t(foo,class_t);
 mp.Digits(50);

% Definiendo Variables simbolicas
% *Nota la parte sumbolica puede ser ejecutada una sola vez*
% centros = x2,y2 
% nodos = x1,y1
sym_vars = {'x1','y1','x2','y2',...
    'NX1','NY1',...
    'NX2','NY2',...
    'dtt','c','muu'};
    
    
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
Sigma = @(f,x,y) (muu)*(jacobian([f(1) f(2)],[x,y]) ...
                    + transpose(jacobian([f(1) f(2)],[x,y])));
                
Sigman = @(f,x,y,nX,nY)  Sigma(f,x,y)*[nX ; nY];
SlipB = @(f,x,y,nX,nY)  [-nY,nX]*Sigman(f,x,y,nX,nY);
 

L1c = @(f) (f(1)-dtt*muu*laplacian(f(1),[x2,y2]) +dtt*diff(f(3),'x2'));
L2c = @(f) (f(2)-dtt*muu*laplacian(f(2),[x2,y2]) +dtt*diff(f(3),'y2'));
B1c = @(f) f(1)*NX2 + f(2)*NY2;
Bslip_c = @(f) SlipB(f,x2,y2,NX2,NY2);

B1 = @(f) f(1)*NX1 + f(2)*NY1;
Bslip = @(f) SlipB(f,x1,y1,NX1,NY1);
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
OperadoresY = {L1c,L2c,B1c,Bslip_c}; %operadores para el anzatz
OperadoresX = {L1,L2,B1,Bslip,Px,Py}; %operadores que le vamos aplicar al anzatz

% Realizando el calculo simbolico
disp('comenzando el calculo symbolico')
symbolic_compute_wolfram(SaveDire,sym_vars,fbr,OperadoresY,OperadoresX);
disp('fin del calculo symbolico')
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
NOperadoresY = length(OperadoresY);
NoperadoresX = length(OperadoresX);
NewDir = ['ConstR',SaveDire];
mkdir(NewDir);
for i=1:NOperadoresY
    for j=1:NoperadoresX
        GenericName = 'fbrGram';
        name = [GenericName,num2str(i),num2str(j)];
        FBRfile = fileread([SaveDire,name,'.m']);
        Patron_de_constantes  = '([0-9]*[.])+[0-9]+([e][+-][0-9])?';
        matchc = regexp(FBRfile,Patron_de_constantes,'match');
        Constantes = unique(matchc);
        
        FBRfilenew = FBRfile;
        for l =1:length(Constantes)
            ConstMp = ['mp(',char(39),Constantes{l},char(39),')'];
            ExpresionRepla = ['(?<!\d)(',replace(Constantes{l},{'.','+','e'},{'[.]','[+]','[e]'}),')(?!(\d+[e]?|[e]))'];
            FBRfilenew = regexprep(FBRfilenew,ExpresionRepla,ConstMp);
        end
        
        fileID = fopen([NewDir,name,'.m'],'w');
        fprintf(fileID,'%s',FBRfilenew);
        fclose(fileID);
    end
end

for i=1:3
    for j=1:NOperadoresY
        GenericName = 'fbrAnzatz';
        name = [GenericName,num2str(i),num2str(j)];
        FBRfile = fileread([SaveDire,name,'.m']);
        Patron_de_constantes  = '([0-9]*[.])+[0-9]+([e][+-][0-9])?';
        matchc = regexp(FBRfile,Patron_de_constantes,'match');
        Constantes = unique(matchc);
        
        FBRfilenew = FBRfile;
        for l =1:length(Constantes)
            ConstMp = ['mp(',char(39),Constantes{l},char(39),')'];
            ExpresionRepla = ['(?<!\d)(',replace(Constantes{l},{'.','+','e'},{'[.]','[+]','[e]'}),')(?!(\d+[e]?|[e]))'];
            FBRfilenew = regexprep(FBRfilenew,ExpresionRepla,ConstMp);
        end
        
        fileID = fopen([NewDir,name,'.m'],'w');
        fprintf(fileID,'%s',FBRfilenew);
        fclose(fileID);
    end
end

addpath('locallib/')
Symbolic_normal_angular(SaveDire)
nombres={'AngNormalx.m','AngNormaly.m'};
    for j=1:length(nombres)
        FBRfile = fileread([SaveDire,nombres{j}]);
        Patron_de_constantes  = '([0-9]*[.])+[0-9]+([e][+-][0-9])?';
        matchc = regexp(FBRfile,Patron_de_constantes,'match');
        Constantes = unique(matchc);
        
        FBRfilenew = FBRfile;
        for l =1:length(Constantes)
            ConstMp = ['mp(',char(39),Constantes{l},char(39),')'];
            ExpresionRepla = ['(?<!\d)(',replace(Constantes{l},{'.','+','e'},{'[.]','[+]','[e]'}),')(?!(\d+[e]?|[e]))'];
            FBRfilenew = regexprep(FBRfilenew,ExpresionRepla,ConstMp);
        end
        
        fileID = fopen([NewDir,nombres{j}],'w');
        fprintf(fileID,'%s',FBRfilenew);
        fclose(fileID);
    end

