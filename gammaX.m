% gammaX.m
% multiply gammaX by input vector/matrix
% 
%           (0  0  0  1)
% gamma1 =  (0  0  1  0) 
%           (0  1  0  0) 
%           (1  0  0  0) 
%
%           (0   0  0  -i)
% gamma2 =  (0   0  i   0) 
%           (0  -i  0   0) 
%           (i   0  0   0) 
%
%           (0  0   1   0)
% gamma3 =  (0  0   0  -1) 
%           (1  0   0   0) 
%           (0  -1  0   0) 
%
%           (1  0  0   0)
% gamma4 =  (0  1  0   0) 
%           (0  0  -1  0) 
%           (0  0  0  -1) 
%
%           (0  0  -i  0)
% gamma5 =  (0  0  0  -i) 
%           (i  0  0   0) 
%           (0  i  0   0) 
%
%
% Assumptions: number of colors = 3
%              number of dirac indicies = 4
%              positions=nx*ny*nz*nt
%           therefore rank =nx*ny*nz*nt*nc*nd
% Input: X   = which gamma to use, X={1,2,3,4,5}
%        inM = input matrix or vector which is [rank x #cols]

function [outM] = gammaX(inM, X)

outM=zeros(size(inM));
n    = size(inM,1);
cols = size(inM,2);

switch X
  case 1  % gamma1
    outM(1:4:n,1:cols) = inM(4:4:n,1:cols);
    outM(2:4:n,1:cols) = inM(3:4:n,1:cols);
    outM(3:4:n,1:cols) = inM(2:4:n,1:cols);
    outM(4:4:n,1:cols) = inM(1:4:n,1:cols);
  case 2  % gamma2
    outM(1:4:n,1:cols) = -1i*inM(4:4:n,1:cols);
    outM(2:4:n,1:cols) =  1i*inM(3:4:n,1:cols);
    outM(3:4:n,1:cols) = -1i*inM(2:4:n,1:cols);
    outM(4:4:n,1:cols) =  1i*inM(1:4:n,1:cols);
  case 3  % gamma3
    outM(1:4:n,1:cols) =  inM(3:4:n,1:cols);
    outM(2:4:n,1:cols) = -inM(4:4:n,1:cols);
    outM(3:4:n,1:cols) =  inM(1:4:n,1:cols);
    outM(4:4:n,1:cols) = -inM(2:4:n,1:cols);
  case 4  % gamma4
    outM(1:4:n,1:cols) =  inM(1:4:n,1:cols);
    outM(2:4:n,1:cols) =  inM(2:4:n,1:cols);
    outM(3:4:n,1:cols) = -inM(3:4:n,1:cols);
    outM(4:4:n,1:cols) = -inM(4:4:n,1:cols);
  case 5  % gamma5
    outM(1:4:n,1:cols) = -1i*inM(3:4:n,1:cols);
    outM(2:4:n,1:cols) = -1i*inM(4:4:n,1:cols);
    outM(3:4:n,1:cols) =  1i*inM(1:4:n,1:cols);
    outM(4:4:n,1:cols) =  1i*inM(2:4:n,1:cols);
  otherwise
    outM = 0;
end
