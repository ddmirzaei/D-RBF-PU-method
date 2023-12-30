function P = PolyMat(X,op)
% This function computes the polynomial matrix (Vandermonde) on a set point
% Inputs:
%   X: points coordinates of size (N x d)
%   op: operator (x,y,xx,yy,L,...)
% Outputs:
%   P: Vandermonde matrix of size (N x Q) on monomials at X

% Note: This function works in 2D for polynomial orders <= 5,
%   to generalize to 3D and higher orders just change the MultiIndex vector

global RBFinfo
PolyOrder = RBFinfo.poly; 
dim = size(X,2);
if PolyOrder < 1
    P = [];
elseif PolyOrder > 5 || dim > 2
    error('for polynomial orders > 5 and dimension > 2 extend the MultiIndex vector in "PolyMat" function')
else
StrVec = @(j,d) [zeros(j-1,1);1;zeros(d-j,1)];
MultiIndex= [0 0;1 0;0 1;2 0;1 1;0 2;3 0;2 1;1 2;0 3;4 0;3 1;2 2;1 3;0 4];
Q = nchoosek(PolyOrder-1+dim,dim); % dimension of polynomial space in R^d
P = [];
switch op
    case '1'
        for k = 1:Q
            P = [P prod(X.^MultiIndex(k,:),2)];
        end
    case 'x'
        for k = 1:Q
            A = MultiIndex(k,:)-StrVec(1,dim)';
            a1 = MultiIndex(k,1); Pc = a1*prod((X+eps).^A,2);
            P = [P Pc];
        end
    case 'y'
        for k = 1:Q
            A = MultiIndex(k,:)-StrVec(2,dim)';
            a2 = MultiIndex(k,2); Pc = a2*prod((X+eps).^A,2);
            P = [P Pc];
        end
    case 'xx'
        for k = 1:Q
            A = MultiIndex(k,:)-2*StrVec(1,dim)';
            a1 = MultiIndex(k,1); Pc = a1*(a1-1)*prod((X+eps).^A,2);
            P = [P Pc];
        end
    case 'yy'
        for k = 1:Q
            A = MultiIndex(k,:)-2*StrVec(2,dim)';
            a2 = MultiIndex(k,2); Pc = a2*(a2-1)*prod((X+eps).^A,2);
            P = [P Pc];
        end
    case 'xy'
        for k = 1:Q
            A = MultiIndex(k,:)-StrVec(1,dim)'-StrVec(2,dim)';
            a1 = MultiIndex(k,1); a2 = MultiIndex(k,2); 
            Pc = a1*a2*prod((X+eps).^A,2);
            P = [P Pc];
        end        
    case 'L' % laplacian
        Pxx = PolyMat(X,'xx'); Pyy = PolyMat(X,'yy');
        P = Pxx + Pyy;
    otherwise
        error('This type of derivative operator op is not implemented in PolyMat function')
end
end