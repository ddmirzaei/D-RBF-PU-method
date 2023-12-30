function w = PUweight(X,Y,rho)
% This function uses to compute the Shepard weights at X where patches 
%   with centers Xc contributes  
% Inouts:
%   X: evaluation points
%   Xc: contributed patches
%   rho: patch radii vector
% Outputs:
%   w: weight vector
%

n = size(X,1);
r = DistMat(X,Y);
if length(rho) == 1
    rho = rho*ones(size(Y,1),1);
end
d = repmat(rho',n,1);
r = r./d;
w = max(1-r,0).^4.*(4*r + 1);
