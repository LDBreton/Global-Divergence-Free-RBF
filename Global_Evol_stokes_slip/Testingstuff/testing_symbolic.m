syms x1 y1 x3


g = x1 + y1 + x3;
matlabFunction(g,'File','testinggg','vars',[x1,y1]);

Symbolic_normal_angular('NormalAngular/')
testinggg(1,1)