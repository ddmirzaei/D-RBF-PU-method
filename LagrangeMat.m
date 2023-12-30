function Lmat = LagrangeMat(Xe,X,xc,rho,op,ScaleOrd)

global RBFinfo 

RBFscale  = RBFinfo.scale;
do_scaling = RBFinfo.do_scaling;

% This function computes the RBF Lagrange functions  
%  for polyharmonic splines (PHS) kernels by applying the scaling rule 
% Inputs:
%   Xe: evaluation points of size ne 
%   X: trial points (centers) of size nx
%   xc: the center of stencil
%   rho: the size of stencil
%   op: operator
%   ScaleOrd: the scaling order of op
% Outputs:
%   Lmat: Lagrange matrix at Xe bases on X

numop = length(op);
Lmat = cell(1,numop);
if do_scaling
    X = (X-xc)/rho;
    Xe= (Xe-xc)/rho;
else
    X = X-xc;
    Xe = Xe-xc;
    ScaleOrd = zeors(1,numop);
end

P = PolyMat(X,'1');
np = size(P,2);
A = KerMat(X,X, '1');
AP = [A P;P' zeros(np,np)];
[LL,UU,PP] = lu(AP);
for k=1:numop
    R = [KerMat(Xe, X, op{k})'; PolyMat(Xe,op{k})'];
    opts.LT = true;
    yy = linsolve(LL,PP*R,opts);
    opt.UT = true;
    Lm = linsolve(UU,yy,opt);
    Lmat{k} = Lm(1:end-np,:)'/rho^ScaleOrd(k);
end

