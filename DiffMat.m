function DM = DiffMat(Xe, X)
% This function calculates the matrix of point differences
% Inputs:
%  Xe: evaluation points of size (m x dim)
%  X: trial points of size (n x dim)
% Outputs:
%  DM: difference matrix of size (m x n)

[dr,cc]=ndgrid(Xe(:),X(:));
DM=dr-cc;