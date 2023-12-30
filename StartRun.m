%%   
% We consider an equation of the form (Lu = -Delta u + u = f) in 2D
% on domain [0,1]^2 with Dirichlet BC (left and right) and Neumann BC (top and bottom)
% True solution: Franke's function 
%           u_y = g2
%         ------------
%        |            |
%        |            |
% u = g1 |   Lu = f   | u = g1
%        |            |
%        |            |
%         ------------
%          u_y = g2

%%
clc
clear

global RBFinfo PUweightType

PointType = 'halton';         % halton or grid                               
RBFtype = 'tp';             % tp (TPS r^k log r, k even) or p (power r^k, k odd)
PUweightType = 'ConstGen';  % Smooth or ConstGen
RBFinfo.type = RBFtype;
RBFinfo.scale = 1;      % not improtant for PHS kernels tp and p
RBFinfo.do_scaling = 1; % 1 for PHS kernels and 0 for other kernels 


disp('-------------------------------------------------------')
disp('D-RBF-PU method strats ...')
if strcmp(RBFinfo.type, 'tp')
fprintf('PUweithType = %7s, RBFtype = r^k log r, (k = 4,6,8)\n',PUweightType) 
else
fprintf('PUweithType = %7s, RBFtype = r^k, (k = 5,7,9)\n',PUweightType) 
end


nlevel = 8;

for p = 1:3
    if strcmp(RBFinfo.type, 'tp')
       RBFinfo.par = (p+1)*2 ; RBFinfo.poly = RBFinfo.par/2+1; % exponents 4,6,8
    else
       RBFinfo.par = (p+1)*2+1; RBFinfo.poly = ceil(RBFinfo.par/2); % exponents 5,7,9 
    end
    h = 0.1;
    for k = 1:nlevel        
        % X: trial points
        % XI: internal test points
        % XB: all boundary points
        % Xb: boundary points on the bottom side of the square
        % Xt: boundary points on the top side of the square
        % Xl: boundary points on the left side of the square
        % Xr: boundary points on the right side of the square
        hh(k) = h;
        [X,XI,XB,Xb,Xt,Xl,Xr] = ScatPoints2D(0,1,0,1,h,PointType);       
        N(k)=size(X,1); NI = size(XI,1); Nb = size(Xb,1); Nt = size(Xt,1);
        
        % patch centers (grid points):
        C_cov = 3;
        hc = C_cov*h; n_min = ceil(0.8*pi*C_cov^2);
        nc = ceil(sqrt(N(k)/C_cov^2)); x = linspace(C_cov/2*h,1-C_cov/2*h,nc);
        [x,y] = meshgrid(x,x); 
        Xc = [x(:) y(:)];
        
        tic
        A = D_RBF_PU(XI,X,Xc,n_min,{'L','1'}); 
        AL = A{1}; A1 = A{2};
        By = D_RBF_PU([Xb; Xt],X,Xc,n_min,'y');
        B1 = D_RBF_PU([Xl;Xr],X,Xc,n_min,'1');
        K = [-AL + A1; By; B1];
        f = -ExactFunc(XI,'L')+ ExactFunc(XI,'1');
        g1 = ExactFunc([Xb; Xt],'y');
        g2 = ExactFunc([Xl; Xr],'1');
        rhs = [f; g1; g2];

        SetupTime(k,p) = toc;  % setup time

        tic
        Uap = K\rhs; 
        SolveTime(k,p) = toc;  % solving time

        InfErr(k,p) = norm(Uap-ExactFunc(X,'1'),inf)./norm(ExactFunc(X,'1'),inf);
        nz_percent(k,p) = nnz(K)/prod(size(K))*100;
        
        h = h/sqrt(2);  % point refinement         
    end    
    % computational orders using a linear curve fitting to errors
    a = polyfit(log(sqrt(N')),log(InfErr(:,p)),1);
    Orders(p) = -a(1);
end
nz_percent
InfErr
Orders

figure;
loglog(sqrt(N),InfErr(:,1),'rs-','MarkerSize',6,'MarkerFaceColor','r','LineWidth',1.3)
hold on
loglog(sqrt(N),InfErr(:,2),'bo-','MarkerSize',6,'MarkerFaceColor','b','LineWidth',1.3)
hold on
loglog(sqrt(N),InfErr(:,3),'m>-','MarkerSize',6,'MarkerFaceColor','m','LineWidth',1.3)
xlabel('$\sqrt{N}$','interpreter','latex')
ylabel('$\|e\|_\infty/\|u\|_\infty$','interpreter','latex')
set(gca,'TickLabelInterpreter','latex')
leg = legend('PHS+P2','PHS+P3', 'PHS+P4','Location','northeast');
set(leg,'Interpreter','latex');
set(gcf, 'Position', [300 300 350 400])
title('D-RBF-PU','Interpreter','latex');
hold off

