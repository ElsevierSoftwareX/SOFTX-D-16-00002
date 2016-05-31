% Author: Mariette Annergren & Christian A. Larsson
% Copyright (c) 2015 Mariette Annergren & Christian A. Larsson

function [V,y] = Vapp(theta,theta0,y0,r,Ny,Nu,Ub,Yb)
   x = 0;
   u = 0;
   y = zeros(size(y0));
   for tk = 1:numel(y)
      y(tk) = theta0(1)*x;
      u     = simpleMPC(theta(2),1,theta(1),Ny,Nu,0,x,u,r(tk),Ub,Yb);
      x     = theta0(2)*x + u;
   end
   V = 1/numel(y0)*(y-y0)'*(y-y0);
end