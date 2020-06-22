function [chi_s]    =   volscat2(tts,ttli)

%Volscatt version 2.
%created by W. Verhoef
%edited by Joris Timmermans to matlab nomenclature.
% date: 11 February 2008
%tts    [1]         Sun            zenith angle in degrees
%tto    [1]         Observation    zenith angle in degrees
%psi    [1]         Difference of  azimuth angle between solar and viewing position
%ttli   [ttli]      leaf inclination array

deg2rad = pi/180;
% nli     = length(ttli);

cos_ttli        = cos(ttli*deg2rad);                %   cosine of normal of upperside of leaf
sin_ttli        = sin(ttli*deg2rad);                %   sine   of normal of upperside of leaf

cos_tts         = cos(tts*deg2rad);                 %   cosine of sun zenith angle
sin_tts         = sin(tts*deg2rad);                 %   sine   of sun zenith angle

Cs              = cos_ttli*cos_tts;                 %   p305{1}
Ss              = sin_ttli*sin_tts;                 %   p305{1}                          


As              = max([Ss,Cs],[],2);

bts             = acos(-Cs./As);                    %   p305{1}

chi_s           = 2/pi*((bts-pi/2).*Cs + sin(bts).*Ss);
