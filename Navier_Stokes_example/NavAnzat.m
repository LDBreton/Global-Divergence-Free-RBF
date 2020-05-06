function [fbrGram2,fbrAnzatz2] = NavAnzat(fbrGram,fbrAnzatz)

[N,M] = size(fbrGram);

fbrGram2 = cell(size(fbrGram));
for i = 1:N
    for j = 1:M
    aux = @(x1,y1,x2,y2,U1,U2) fbrGram{i,j}(x1,y1,x2,y2,U1(double(x2),double(y2)),U2(double(x2),double(y2)),...
                                                       U1(double(x1),double(y1)),U2(double(x1),double(y1)) );
    fbrGram2(i,j) = {aux};
    end
end

[N,M] = size(fbrAnzatz);

fbrAnzatz2 = cell(size(fbrAnzatz));
for i = 1:N
    for j = 1:M
    aux = @(x1,y1,x2,y2,U1,U2) fbrAnzatz{i,j}(x1,y1,x2,y2,U1(double(x2),double(y2)),U2(double(x2),double(y2)),...
                                                       U1(double(x1),double(y1)),U2(double(x1),double(y1)) );
    fbrAnzatz2(i,j) = {aux};
    end
end

end

