muu=1.0e-06;
FileName=['Tabletesting_wolfram/R2Errores_muu_',num2str(muu),'_dirichbdf2.txt'];
Dirsave = '/home/porufes/Dropbox/version_editorial-5-mayo/version_octubre/tables';
Fileout = [Dirsave,'/Table_dirich_muu',num2str(muu),'.tex'];
caption = ['Dirichlet Boundary condition with $\mu$ = ',num2str(muu)];

Data = readmatrix(FileName);
VecTimes = unique(Data(:,2));
VecTimes = VecTimes(end:-1:1);
VecNodos = unique(Data(:,1));

Ntime = length(VecTimes);
ErroresY = reshape(Data(:,3),[],Ntime);
ErroresP = reshape(Data(:,4),[],Ntime);
CondA = reshape(Data(:,5),Ntime,[])';

LatexTables = zeros(size([ErroresY,ErroresP]));

for i=1:length(VecNodos)
LatexTables(i,:) = reshape([ErroresY(i,:);ErroresP(i,:)],1,[]);
end


fileID = fopen(Fileout,'w');
% fprintf(fileID,'\\begin{table} \n \\begin{centering} \n ');
% fprintf(fileID,'\\begin{tabular}{%s',repmat('l',1,size(LatexTables,2)+2));
fprintf(fileID,'\\begin{tabular}{l|');
fprintf(fileID,'%s',repmat('ll|',1,Ntime));
fprintf(fileID,'l} \n');
fprintf(fileID,'$N\\backslash\\Delta t$');
fprintf(fileID,' & \\multicolumn{2}{c}{%.2e}',VecTimes);
fprintf(fileID,' & Max Cond \\tabularnewline \n');
labels = '& $||e_{y}||_{\infty}$ & $||e_{\nabla p}||_{\infty}$ ';
fprintf(fileID,'%s',repmat(labels,1,Ntime));
fprintf(fileID,' & \\tabularnewline \n');

for i=1:length(VecNodos)
fprintf(fileID,' %d ',VecNodos(i));
fprintf(fileID,' & %.2e',LatexTables(i,:));
fprintf(fileID,'& %.2e \\tabularnewline \n',max(CondA(i,:)));
end
fprintf(fileID,'\\end{tabular} \n');
% fprintf(fileID,'\\par\\end{centering} \n');
% fprintf(fileID,'\\caption{%s} \n',caption);
% fprintf(fileID,'\\end{table}');

fclose(fileID);