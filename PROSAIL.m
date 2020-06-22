function [Gs,Go,k,K,sob,sof]    =   PROSAIL(tts,tto,psi,lidf)

deg2rad = pi/180;

litab=[5.,15.,25.,35.,45.,55.,65.,75.,81.,83.,85.,87.,89.]';
%% 1.0 Geometric quantities
% 1.1 general geometric quantities
% these variables are scalars
cos_tts     = cos(tts*deg2rad);             %           cos solar       angle   
% tan_tto     = tan(tto*deg2rad);             %           tan observation angle

cos_tto     = cos(tto*deg2rad);             %           cos observation angle   
% sin_tts     = sin(tts*deg2rad);             %           sin solar       angle

% tan_tts     = tan(tts*deg2rad);             %           tan solar angle

psi         = abs(psi-360*round(psi/360));  %           (to ensure that volscatt is symmetric for psi=90 and psi=270)
% dso         = sqrt(tan_tts.^2 + tan_tto.^2 - 2*tan_tts.*tan_tto.*cos(psi*deg2rad));
% 1.2 geometric factors associated with extinction and scattering
[chi_s,chi_o,frho,ftau]=volscat(tts,tto,psi,litab);   % volume scattering

ksli        = chi_s./cos_tts;               % [13]      p306{1} extinction coefficient in direction of sun        per leaf angle
koli        = chi_o./cos_tto;               % [13]      p307{1} extinction coefficient in direction of observer   per leaf angle

sobli       = frho*pi/(cos_tts*cos_tto);    % [13]      pag 309{1} area scattering coefficient fractions
sofli       = ftau*pi/(cos_tts*cos_tto);    % [13]      pag 309{1}

%integration over angles (using a vector inproduct) -> scalars
k           = ksli'*lidf;                   %           pag 306{1}    extinction coefficient in direction of sun.
K           = koli'*lidf;                   %           pag 307{1}    extinction coefficient in direction of observer
sob         = sobli'*lidf;                  %           weight of specular2directional back    scatter coefficient
sof         = sofli'*lidf;                  %           weight of specular2directional forward scatter coefficient
% wso           = sob*rho + sof*tau;            % [nwl]     w,      p309{1} bidirectional scattering coefficent (directional-directional)  

% 1.5 probabilities Ps, Po, Pso
% Ps          =   exp(-k*LAI);                                              % [nl+1]  p154{1} probability of viewing a leaf in solar dir
% Po          =   exp(-K*LAI);                                              % [nl+1]  p154{1} probability of viewing a leaf in observation dir

Gs=k*cos_tts;
Go=K*cos_tto;

% i0=1-Ps;

