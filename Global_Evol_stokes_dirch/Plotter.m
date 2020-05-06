muu=1e-06;
FileName=['Errores_muu_',num2str(muu),'_dirich.txt'];
Dirsave = '/home/porufes/Dropbox/version_editorial-5-mayo/version_octubre/Image';
Fileout = [Dirsave,'/Dirich_muu_y',num2str(muu)];
indentifier = 'Dirichlet';
k=3;
LimY = 0.3;



if(k==4)
YLabel = '$||e_{\nabla p}||_{\infty}$ ';
else
YLabel = '$||e_{y}||_{\infty}$ ';
end    

Data = readmatrix(FileName);
VecTimes = unique(Data(:,2));
VecTimes = VecTimes(end:-1:1);
VecNodos = unique(Data(:,1));
Ntime = length(VecTimes);
% ErroresY = reshape(Data(:,3),Ntime,[])';
% ErroresP = reshape(Data(:,4),Ntime,[])';

Dataplot = reshape(Data(:,k),Ntime,[])';

titulo = [YLabel,indentifier,' $\mu$ =',num2str(muu)];
Markers = {'+','*','x','v','d','^','s','>','<'};
Linest ={'-','--',':','-.'};
figure1 = figure(1);
hold on
Leng = {};
for i=2:length(VecNodos)
loglog(VecTimes, Dataplot(i,:), [Linest{mod(i,4)+1},Markers{i}],'MarkerSize',10,'LineWidth',3);
Leng = [Leng,{['N ',num2str(VecNodos(i))]}];
end
hold off
xlabel('time step');
ylabel(YLabel,'Interpreter','latex');
xtickformat('%0.1e')
ytickformat('%1.3f')
set(gca,'XTickLabelRotation',45)
ylim([0,LimY])
xticks('auto')
grid on
set(gca,'FontSize',30)
legend(Leng)
title(titulo,'Interpreter','latex')
figure1.Position=[66 1 1615 971];
saveas(figure1,Fileout,'epsc');
savefig(Fileout)
