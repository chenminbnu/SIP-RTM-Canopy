function [chi_s,chi_o,frho,ftau]    =   volscat(tts,tto,psi,ttli)

%Volscatt version 2.
%created by W. Verhoef
%edited by Joris Timmermans to matlab nomenclature.
% date: 11 February 2008
%tts    [1]         Sun            zenith angle in degrees
%tto    [1]         Observation    zenith angle in degrees
%psi    [1]         Difference of  azimuth angle between solar and viewing position
%ttli   [ttli]      leaf inclination array

deg2rad = pi/180;
nli     = length(ttli);

psi_rad         = psi*deg2rad*ones(nli,1);

cos_psi         = cos(psi*deg2rad);                 %   cosine of relative azimuth angle

cos_ttli        = cos(ttli*deg2rad);                %   cosine of normal of upperside of leaf
sin_ttli        = sin(ttli*deg2rad);                %   sine   of normal of upperside of leaf

cos_tts         = cos(tts*deg2rad);                 %   cosine of sun zenith angle
sin_tts         = sin(tts*deg2rad);                 %   sine   of sun zenith angle

cos_tto         = cos(tto*deg2rad);                 %   cosine of observer zenith angle
sin_tto         = sin(tto*deg2rad);                 %   sine   of observer zenith angle

Cs              = cos_ttli*cos_tts;                 %   p305{1}
Ss              = sin_ttli*sin_tts;                 %   p305{1}                          

Co              = cos_ttli*cos_tto;                 %   p305{1}
So              = sin_ttli*sin_tto;                 %   p305{1}

As              = max([Ss,Cs],[],2);
Ao              = max([So,Co],[],2);

bts             = acos(-Cs./As);                    %   p305{1}
bto             = acos(-Co./Ao);                    %   p305{2}

chi_o           = 2/pi*((bto-pi/2).*Co + sin(bto).*So);
chi_s           = 2/pi*((bts-pi/2).*Cs + sin(bts).*Ss);

delta1          = abs(bts-bto);                     %   p308{1}
delta2          = pi-abs(bts + bto - pi);           %   p308{1}

Tot             = psi_rad + delta1 + delta2;        %   pag 130{1}

bt1             = min([psi_rad,delta1],[],2);
bt3             = max([psi_rad,delta2],[],2);
bt2             = Tot - bt1 - bt3;

T1              = 2.*Cs.*Co + Ss.*So.*cos_psi;
T2              = sin(bt2).*(2*As.*Ao + Ss.*So.*cos(bt1).*cos(bt3));

Jmin            = (   bt2).*T1 - T2;
Jplus           = (pi-bt2).*T1 + T2;

frho            =  Jplus/(2*pi^2);
ftau            = -Jmin /(2*pi^2);

% pag.309 wl-> pag 135{1}
frho            = max([zeros(nli,1),frho],[],2);
ftau            = max([zeros(nli,1),ftau],[],2);
