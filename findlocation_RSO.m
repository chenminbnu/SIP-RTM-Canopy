function [location]=findlocation_RSO(INPUT)

i_tts= min(max(round((INPUT(:,1)-5)/10)+1,1),7);
i_tto= min(max(round((INPUT(:,2))/10)+1,1),7);
INPUT(:,3)=abs(INPUT(:,3)-360*round(INPUT(:,3)/360)); 
i_psi= min(max(round((INPUT(:,3))/30)+1,1),7);
i_CIy1=round((INPUT(:,4)-0.4)*10);
i_CIy2=round((INPUT(:,5)-0.4)*10);
i_LAI=min(max(round(INPUT(:,6)*10),1),70);
i_LIDFa=INPUT(:,7);
i_rsoil=INPUT(:,8);
location=(i_tts-1)*7*7*6*6*70*8*10+(i_tto-1)*7*6*6*70*8*10+(i_psi-1)*6*6*70*8*10+(i_CIy1-1)*6*70*8*10+(i_CIy2-1)*70*8*10+(i_LAI-1)*8*10+(i_LIDFa-1)*10+i_rsoil;
end
