function Ind = PointsInPatch(X,Xc,r)
% This function collects the indices of points in X in patches in Xc
%  using the kd-tree algorithm 
% Inputs:
%   X: a set of points 
%   Xc: patch's centers
%   r : radius of patches
% Outputs:
%   Ind: indices of points in each patch

T = KDTreeSearcher(X);   
Ind = rangesearch(T,Xc,r);
