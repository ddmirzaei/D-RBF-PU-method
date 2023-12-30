function ScaleOrd = ScalingOrder (op)
% This function returns the scaling order of operator op
% Inputs:
%  op: operator
% Outputs:
%  ScaleOrd: the scaling order
%
numop = length(op);
for k = 1:numop
switch op{k}
    case '1', s = 0; case 'x', s = 1; case 'y', s = 1;
    case 'xx', s = 2;case 'xy',s = 2; case 'yy',s = 2;
    case 'L', s = 2;
    otherwise
        error('The scaling order of operator op is not implemented')
end
ScaleOrd(k) = s;
end
