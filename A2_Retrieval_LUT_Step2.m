% 
% load RSO_LUT.mat 
% load RDO_LUT.mat 
% 
% load RSD_LUT.mat 
% load RDD_LUT.mat 


%Directional view
RO_LUT=zeros(11*7*7*7*6*6*70*8,8);
for i_Ratio=1:11    
    k2=(i_Ratio-1)*0.1;
    k1=1-k2;    
    for i_tts=1:7
        for i_tto=1:7
%             for i_psi=1:7
%                 for i_CIy1=1:6
%                     for i_CIy2=1:6
%                         for i_LAI=1:70
%                             for i_LIDFa=1:8                             

                                RO_LUT((i_Ratio-1)*7*7*7*6*6*70*8+(i_tts-1)*7*7*6*6*70*8+(i_tto-1)*7*6*6*70*8+1:(i_Ratio-1)*7*7*7*6*6*70*8+(i_tts-1)*7*7*6*6*70*8+i_tto*7*6*6*70*8,:)=k1*RSO_LUT((i_tts-1)*7*7*6*6*70*8+(i_tto-1)*7*6*6*70*8+1:(i_tts-1)*7*7*6*6*70*8+i_tto*7*6*6*70*8,:)+k2*repmat(RDO_LUT((i_tto-1)*6*6*70*8+1:i_tto*6*6*70*8,:),7,1);
                                
%                             end
%                         end
%                     end
%                 end
%             end
        end
    end        
end
RO_LUT(:,8)=(RO_LUT(:,2)-RO_LUT(:,1))./(RO_LUT(:,2)+RO_LUT(:,1)).*RO_LUT(:,2);
save RO_LUT.mat RO_LUT;



%Hemi view
RD_LUT=zeros(11*7*6*6*70*8,8);
for i_Ratio=1:11    
    k2=(i_Ratio-1)*0.1;
    k1=1-k2;    
%     for i_tts=1:7
%                 for i_CIy1=1:6
%                     for i_CIy2=1:6
%                         for i_LAI=1:70
%                             for i_LIDFa=1:8                             

                                RD_LUT((i_Ratio-1)*7*6*6*70*8+1:i_Ratio*7*6*6*70*8,:)=k1*RSD_LUT+k2*repmat(RDD_LUT,7,1);
                                
%                             end
%                         end
%                     end
%                 end
%     end        
end
RD_LUT(:,8)=(RD_LUT(:,2)-RD_LUT(:,1))./(RD_LUT(:,2)+RD_LUT(:,1)).*RD_LUT(:,2);
save RD_LUT.mat RD_LUT;
    
    