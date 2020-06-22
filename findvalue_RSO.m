function [OUTPUT]=findvalue_RSO(location)

i_CIy1=floor(location./(6*70*8*10))+1;
i_CIy2=floor((location-(i_CIy1-1)*6*70*8*10)./(70*8*10))+1;
i_LAI=floor((location-(i_CIy1-1)*6*70*8*10-(i_CIy2-1)*70*8*10)./(8*10))+1;
i_LIDFa=floor((location-(i_CIy1-1)*6*70*8*10-(i_CIy2-1)*70*8*10-(i_LAI-1)*8*10)./10)+1;
i_Rsoil =location-(i_CIy1-1)*6*70*8*10-(i_CIy2-1)*70*8*10-(i_LAI-1)*8*10-(i_LIDFa-1)*10;
CIy1=(i_CIy1-1)*0.1+0.5;
CIy2=(i_CIy2-1)*0.1+0.5;
LAI=i_LAI*0.1;
LIDFa=i_LIDFa*10;
Rsoil = i_Rsoil;
OUTPUT=[CIy1 CIy2 LAI LIDFa i_Rsoil];
end
