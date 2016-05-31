% Example application cost where the estimate
% of theta(2) is important

% Author: Mariette Annergren & Christian A. Larsson
% Copyright (c) 2015 Mariette Annergren & Christian A. Larsson

function V = Vapp(theta)
    theta0 = [10 -9];
    V      = norm(theta(2)-theta0(2),2);
end