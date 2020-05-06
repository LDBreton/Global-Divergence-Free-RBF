function [fbrGram2] = RemplaceNormal(fbrGram,Nx,Ny)

[N,M] = size(fbrGram);

fbrGram2 = cell(size(fbrGram));
for i = 1:N
    for j = 1:M
    aux = @(x1,y1,x2,y2,varargin) fbrGram{i,j}(x1,y1,x2,y2...
                                   ,Nx(x1,y1)...
                                   ,Ny(x1,y1)...
                                   ,Nx(x2,y2)...
                                   ,Ny(x2,y2)...
                                   ,varargin{:});
    fbrGram2(i,j) = {aux};
    end
end


end