function [CIs]=CIxy(CIy1,CIy2,tts)
if tts<20
    CIs=CIy1;
else if tts>60
    CIs=CIy2;   
    else
CIs=(CIy2-CIy1)/(60-20)*(tts-20)+CIy1;
    end
end