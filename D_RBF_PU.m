function A = D_RBF_PU(Y,X,Xc,n_min,op)
% This function produces the d-rbf-pu matrix to approximate 
%        (Lu)|Y = A*u|X 
%   where L is a linear operator 
% Inputs:
%   Y: test points
%   X: trial points
%   Xc: patch centers
%   n_min: minimum number of points in each patch (to gaurantee the overlap between patches)
%   op: operator
% Outputs:
%   A: the differentiation matrix 

global PUweightType  % PU weight type: Smooth or Constant-generated
global RBFinfo

if ~iscell(op), op = {op}; end
[N,dim] = size(X); Ne = size(Y,1); Nc = size(Xc,1);  numop = length(op); % sizes ...
for k = 1:numop
    A{k} = spalloc(Ne,N,0);       % initial allocation
end

% Obtaining patch radii involves ensuring that n_loc points are positioned within each patch. 
% The parameter n_min, representing the minimum number of points, ensures that the patches overlap.
q = nchoosek(RBFinfo.poly-1+dim,dim); % dimension of poly space in dim-dimension
n_loc = max(3*q,n_min);   % number of trial points in each patch = 3*q

[IndX,dist] = knnsearch(X,Xc,'k',n_loc); % indices and distances of trial points in each patch
rho = max(dist,[],2);                    % patch radii vector
IndXc = PointsInPatch(Xc,Xc,2*max(rho)); % will be uses to quickly find close patch centers to an specific center 
IndY0 = PointsInPatch(Y,Xc,max(rho));    % will be used to accelerate test point search 


ScaleOrd = ScalingOrder (op);       % scaling order of the operator

switch PUweightType
    case 'Smooth'
        for j=1:Nc
            % xc: patch center number j
            xc = Xc(j,:);        
            % Xj: all trial points in patch j
            ind = IndX(j,:); indc = IndXc{j};
            Xj = X(ind,:); Nind = length(ind);           
            % Yj: all test points in patch j
            inde0 = DistMat(Y(IndY0{j},:),xc) <= rho(j);
            inde = IndY0{j}(inde0); 
            Yj = Y(inde,:);            
            if isempty(inde) == 0                
                % w: PU weights at test points 
                pn = PUweight(Yj,xc,rho(j));
                pd = PUweight(Yj,Xc(indc,:),rho(indc));
                w = pn./sum(pd,2);
                % c: weight matrices (Lagrance functions at Yj)
                c = LagrangeMat(Yj,Xj,xc,rho(j),op,ScaleOrd);
                % inserting weights into global matrices and updating
                for k = 1:numop
                    A{k}(inde,ind) = A{k}(inde,ind) + repmat(w,1,Nind).*c{k};
                end
            end
        end
    case 'ConstGen'
        for j=1:Nc   
            % xc: patch center number j
            xc = Xc(j,:);    
            % Xj: all trial points in patch j
            ind = IndX(j,:); indc = IndXc{j};
            Xj = X(ind,:);        
            % finding indices of test points in which xc is their closest center
            inde0 = DistMat(Y(IndY0{j},:),xc) <= rho(j);
            inde0 = IndY0{j}(inde0); 
            Yj0 = Y(inde0,:);
            D = DistMat(Yj0,Xc(indc,:));
            [Dmin,Dind] = min(D,[],2); % find minimum distances           
            inde = (Dind == 1);
            indx = inde0(inde); % going from local to global indices
            if ~isempty(inde)
                % test points in tile j
                Yj = Yj0(inde,:);            
                % c: weight matrices (Lagrance functions at Yj)
                c = LagrangeMat(Yj,Xj,xc,rho(j),op,ScaleOrd);   
                % inserting weights into global matrices
                for k = 1:numop
                    A{k}(indx,ind) = c{k};
                end
            end
        end
end
if numop == 1
    A = A{1};
end

