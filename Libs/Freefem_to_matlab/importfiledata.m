function dat=importfiledata(fileToRead1)
%IMPORTFILE(FILETOREAD1)
%  Imports data from the specified file
%  FILETOREAD1:  file to read

%  Auto-generated by MATLAB on 28-Feb-2010 19:54:44

% Import the file
rawData1 = importdata(fileToRead1);


np=rawData1(1,4);
dat = zeros(1,np);
k=0;
for i=2:np+1
    k=k+1;
    dat(1,k)=rawData1(i,1);
end
