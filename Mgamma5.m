%function to multiply M by gamma5 on the right

function [outM] = Mgamma5(inM,n,cols);

% outM = zeros(n,cols);

outM(1:n,1:4:cols) = 1i*inM(1:n,3:4:cols);
outM(1:n,2:4:cols) = 1i*inM(1:n,4:4:cols);
outM(1:n,3:4:cols) = -1i*inM(1:n,1:4:cols);
outM(1:n,4:4:cols) = -1i*inM(1:n,2:4:cols);