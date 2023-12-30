function [X,XI,XB,Xb,Xt,Xl,Xr] =  PointsInSquare(a,b,c,d,h,ptype)
% This function produces grid or Halton points on square [a,b]x[c,d]
% Inputs:
%   a, b, c, d: the lengths of square sides
%   h: fill distance
%   ptype: type of points, either grid or Halton
% Outputs
    % X: trial points
    % XI: internal test points
    % XB: all boundary points
    % Xb: boundary points on the bottom side of the square
    % Xt: boundary points on the top side of the square
    % Xl: boundary points on the left side of the square
    % Xr: boundary points on the right side of the square
%   
switch ptype
    case 'grid'
        [x1,y1] = meshgrid(a:h:b,c:h:d); X = [x1(:) y1(:)];
        x = a+h:h:b-h; y = c+h:h:d-h;
        [xt,yt] = meshgrid(x,y); XI = [xt(:) yt(:)];
        Xb = [x' c*ones(size(x'))]; Xt = [x' d*ones(size(x'))];
        y = [c y d];
        Xl = [a*ones(size(y')) y']; Xr = [b*ones(size(y')) y'];
        XB = [Xl;Xr;Xb;Xt];
    case 'halton'        
        N = ceil((1/h).^2);
        p = haltonset(2,'Skip',1e4,'Leap',1e2);
        XI = ((1-h)*(2*net(p,N)-1)+1)/2;
        x = a+h:h:b-h; y = c+h:h:d-h;
        Xb =[x' c*ones(size(x'))]; Xt = [x' d*ones(size(x'))];
        y = [c y d];
        Xl = [a*ones(size(y')) y']; Xr = [b*ones(size(y')) y'];
        XB = [Xl;Xr;Xb;Xt];  
        X = [XI;XB];
end
