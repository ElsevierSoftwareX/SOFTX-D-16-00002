% Author: Mariette Annergren & Christian A. Larsson
% Copyright (c) 2015 Mariette Annergren & Christian A. Larsson

function ellipse(varargin)
% FUNCTION:     ellipse
% DESCRIPTION:  Given a matrix A and a point c the function ellipse 
%               draws an ellipse centered in c. [(x+c)' * K * (x+c) = 1]
% INPUT:        1 - matrix [hessian]
%               2 - center point [estimates of true parameters]
%               3 - color 
%               4 - 'empty' [no center point plotted]
% OUTPUT:       plot of ellipse

if nargin == 3
    
    A     = varargin{1}; 
    c     = varargin{2}; 
    a     = varargin{3};
    [Q,D] = eig(A);
    if abs(D(1,1)) <1e-8; D(1,1) = 1e-4; elseif abs(D(2,2)) <1e-8; D(2,2) = 1e-4; end % handles infinite semi-axis
    W     = Q*1/sqrt(D);
    angle = linspace(0,2*pi,100);
    
    hold on
    plot(W(1,1).*cos(angle)+W(1,2).* sin(angle)+c(1),...
         W(2,1).*cos(angle)+W(2,2).* sin(angle)+c(2),...
         a,'LineWidth',2)
    plot(c(1),c(2),'*')    
    axis equal
    
elseif nargin == 4
    
    A     = varargin{1}; 
    c     = varargin{2}; 
    a     = varargin{3};
    [Q,D] = eig(A);
    if abs(D(1,1)) <1e-8; D(1,1) = 1e-4; elseif abs(D(2,2)) <1e-8; D(2,2) = 1e-4; end % handles infinite semi-axis
    W     = Q*1/sqrt(D);
    angle = linspace(0,2*pi,100);
    
    hold on
    plot(W(1,1).*cos(angle)+W(1,2).* sin(angle)+c(1),...
         W(2,1).*cos(angle)+W(2,2).* sin(angle)+c(2),...
         a,'LineWidth',2)
    axis equal
    
end


