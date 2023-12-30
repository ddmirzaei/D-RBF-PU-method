function y =  Frbf (s,k)
% This function computes the k-th derivative of standard PHS kernels
%    in f form, i.e., as a function of s = r^2/2. For more details see
%    "R. Schaback, Matlab programming for kernel based methods, Technical Report, 2011."
%
global RBFinfo;
RBFtype = RBFinfo.type; RBFpar = RBFinfo.par;
switch lower(RBFtype)
    case ('p')            % (-1)^ceil(RBFpar/2)*(r)^(RBFpar)
        beta = RBFpar;
        fac=1;
        for i=0:k-1
            fac=fac*(beta/2-i);
        end
        y = (-1)^ceil(beta/2)*fac*(s+eps).^(beta/2-k);
    case('tp')            % r^RBFpar*log(r)
        ord = k; par=RBFpar;
        fac = 1;
        su = 0;
        while ord>0
            ord=ord-1;
            if (ord== k-1)
                su=1;
            else
                su=fac+su*par;
            end
            fac = fac*par;
            par=par-2;
        end
        y = (2*s+4*eps).^(par/2);
        y=fac*y.*log(2*s+4*eps)/2+su*y;
    otherwise
        error('RBF type not implemented')
end
