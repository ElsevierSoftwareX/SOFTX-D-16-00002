% Application cost emphasizing typical robust control criterion

% Author: Mariette Annergren & Christian A. Larsson
% Copyright (c) 2015 Mariette Annergren & Christian A. Larsson

function Fmax = F(theta,G0,Tmag,wSamp,Ts)

F      = [1 theta(2)];
B      = [0 theta(1)];
G      = tf(B,F,Ts,'Variable','z^-1');      
G_G0   = minreal(G-G0);

Gfr    = squeeze(freqresp(G,wSamp));
G_G0fr = squeeze(freqresp(G_G0,wSamp));

Fmax   = max((G_G0fr'.*G_G0fr.').*Tmag./(Gfr'.*Gfr.'));

end