% gamma5.m
% multiply gamma5 by input vector/matrix
% 
%           (0  0  -1  0)
% gamma5 = i(0  0  0  -1) 
%           (1  0  0   0) 
%           (0  1  0   0) 
%
%
% PRECONDION: positions=nx*ny*nz*nt
%
% Assumptions: number of colors = 3
%              number of dirac indicies = 4

function [outM] = gamma5(inM, n, cols)

outM=zeros(49152,cols);

outM(1:4:n,1:cols) = -1i*inM(3:4:n,1:cols);
outM(2:4:n,1:cols) = -1i*inM(4:4:n,1:cols);
outM(3:4:n,1:cols) =  1i*inM(1:4:n,1:cols);
outM(4:4:n,1:cols) =  1i*inM(2:4:n,1:cols);
