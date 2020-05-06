% Apoint = 1/50;
% Idealproportion = searchprop(Apoint);
% [P_ni,P_f] = Mesh_gen(Apoint,Apoint*Idealproportion,0.5);
%              hold on
%              plot(P_ni(:,1),P_ni(:,2),'.r')
%              plot(P_f(:,1),P_f(:,2),'.b')
%              hold off
% 
% length(P_ni)/length(P_f)


Apoint = 1/25;
[P_ni,P_f] = Mesh_gen(Apoint*0.1,Apoint,0.5);
             hold on
             plot(P_ni(:,1),P_ni(:,2),'.r')
             plot(P_f(:,1),P_f(:,2),'.b')
             hold off