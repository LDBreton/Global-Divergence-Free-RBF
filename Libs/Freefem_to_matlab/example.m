clear;

[points seg tri]=importfilemesh('Gilgamesh.msh');

%dat=importfiledata('Heatb.bb');
dat=importfiledata('Heat.bb');



figure1 = figure(1);
 pdeplot(points,seg,tri,'xydata',dat,'zdata',dat,'mesh','off','colormap','jet');
      view([104 28]);
xlim([min(points(1,:)) max(points(1,:))]);
ylim([min(points(2,:)) max(points(2,:))]);
grid on;    