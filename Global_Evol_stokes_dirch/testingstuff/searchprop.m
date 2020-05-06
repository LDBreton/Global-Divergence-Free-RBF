function Idealproportion = searchprop(Apoint)
func = @(x) prop(Apoint,x);
%Idealproportion = fminbnd(func,0.001,1.0);
options = optimset('Display','iter');
Idealproportion = fzero(func,1.0,options);
end

function Differencefs = prop(Apoint,Idealprop)
[P_ni,P_f] = Mesh_gen(Apoint,Apoint*Idealprop,0.5);
Differencefs = ((length(P_ni)/length(P_f)) - 1);

end