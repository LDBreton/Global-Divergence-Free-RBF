function saveplot(direc,points,seg,tri,j)

dat=importfiledata([direc,'controlsave/controlnorm_',num2str(j),'.txt']);
time = (j+1)*(1/50);
figure1 = figure(1);
pdeplot(points,seg,tri,'xydata',dat,'zdata',dat,'mesh','off','colormap','jet');
view([104 28]);
xlim([min(points(1,:)) max(points(1,:))]);
ylim([min(points(2,:)) max(points(2,:))]);
xlabel('x') % x-axis label
ylabel('y') % y-axis label
title(['||v||_2,time = ', num2str(time)])
grid on;
saveas(figure1,[direc,'control_plot/controlnorm_',num2str(j),'.jpg']);
end