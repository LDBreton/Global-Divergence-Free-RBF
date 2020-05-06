function [P_sc,P_fc] = Mesh_gen_circle_N(N,radius)

command = char(strcat('FreeFem++ -nw -v 0 FreeFem_meshing/Mesh_gen_circle_N.edp', {' '},...
                      '-NN', {' '} ,num2str(N,16),{' '},...
                      '-rad', {' '} ,num2str(radius,16)));
                  
                  
[status,cmdout] = system(command);
Meshdata = textscan(cmdout,'%f %f %d');

idxInsidePoint= (Meshdata{3}==0);

P_sc = [ Meshdata{1}(idxInsidePoint), Meshdata{2}(idxInsidePoint) ];
P_fc = [ Meshdata{1}(~idxInsidePoint), Meshdata{2}(~idxInsidePoint) ];

% hold on
% plot(P_sc(:,1),P_sc(:,2),'.r')
% plot(P_fc(:,1),P_fc(:,2),'.b')
% hold off