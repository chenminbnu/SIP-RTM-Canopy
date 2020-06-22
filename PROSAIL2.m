function [Gs,k]    =   PROSAIL2(tts,lidf)

deg2rad = pi/180;

litab=[5.,15.,25.,35.,45.,55.,65.,75.,81.,83.,85.,87.,89.]';
%% 1.0 Geometric quantities
% 1.1 general geometric quantities
% these variables are scalars
cos_tts     = cos(tts*deg2rad);             %           cos solar       angle   

% sin_tts     = sin(tts*deg2rad);             %           sin solar       angle
% 
% tan_tts     = tan(tts*deg2rad);             %           tan solar angle

% dso         = sqrt(tan_tts.^2 + tan_tto.^2 - 2*tan_tts.*tan_tto.*cos(psi*deg2rad));
% 1.2 geometric factors associated with extinction and scattering
[chi_s]=volscat2(tts,litab);   % volume scattering

ksli        = chi_s./cos_tts;               % [13]      p306{1} extinction coefficient in direction of sun        per leaf angle

%integration over angles (using a vector inproduct) -> scalars
k           = ksli'*lidf;                   %           pag 306{1}    extinction coefficient in direction of sun.

% 1.5 probabilities Ps, Po, Pso
% Ps          =   exp(-k*LAI);                                              % [nl+1]  p154{1} probability of viewing a leaf in solar dir

Gs=k*cos_tts;

% i0=1-Ps;

