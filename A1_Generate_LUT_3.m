% %PART1
% BRFv2_hemi_LUT=zeros(7*6*6*70*8,3);
% 
% ind=0;
% for i_tts=1:7
%     tts=(i_tts-1)*10+5;
%     for i_CIy1=1:6
%         for i_CIy2=1:6
%             CIy1=(i_CIy1-1)*0.1+0.5;
%             CIy2=(i_CIy2-1)*0.1+0.5;
%             [CIs]=CIxy(CIy1,CIy2,tts);
%             for i_LAI=1:70
%                 LAI=i_LAI*0.1;
%                 for i_LIDFa=1:8
%                     LIDFa=i_LIDFa*10;
%                     lidf=campbell(LIDFa);
%                     lidf=lidf';
% 
%                     [sob_vsla, sof_vsla,kgd]=A_BRFv2_single_hemi(tts,CIs,CIy1,CIy2,LAI,lidf);
% 
%                     ind=ind+1;
%                     BRFv2_hemi_LUT(ind,:)=[sob_vsla sof_vsla kgd];
%                 end
%             end
%         end
%     end
% end
% 
% save BRFv2_hemi_LUT.mat BRFv2_hemi_LUT;
% 
% 
% %PART2
% BRFv2_dif_LUT=zeros(7*6*6*70*8,3);
% 
% ind=0;
% for i_tto=1:7
%     tto=(i_tto-1)*10;
%     for i_CIy1=1:6
%         for i_CIy2=1:6
%             CIy1=(i_CIy1-1)*0.1+0.5;
%             CIy2=(i_CIy2-1)*0.1+0.5;
%             [CIo]=CIxy(CIy1,CIy2,tto);
%             for i_LAI=1:70
%                 LAI=i_LAI*0.1;
%                 for i_LIDFa=1:8
%                     LIDFa=i_LIDFa*10;
%                     lidf=campbell(LIDFa);
%                     lidf=lidf';
% 
%                     [sob_vsla_dif, sof_vsla_dif, kg_dif]=A_BRFv2_single_dif(tto,CIo,CIy1,CIy2,LAI,lidf);  %the same for VZA
% 
%                     ind=ind+1;
%                     BRFv2_dif_LUT(ind,:)=[sob_vsla_dif sof_vsla_dif kg_dif];
%                 end
%             end
%         end
%     end
% end
% 
% save BRFv2_dif_LUT.mat BRFv2_dif_LUT;


% PART3
BRFv2_hemi_dif_LUT=zeros(6*6*70*8,3);

ind=0;
for i_CIy1=1:6
    CIy1=(i_CIy1-1)*0.1+0.5;
    for i_CIy2=1:6
        CIy2=(i_CIy2-1)*0.1+0.5;
        for i_LAI=1:70
            LAI=i_LAI*0.1;
            for i_LIDFa=1:8
                LIDFa=i_LIDFa*10;
                lidf=campbell(LIDFa);
                lidf=lidf';
                
                [sob_vsla_hemi_dif, sof_vsla_hemi_dif,kgd_dif]=A_BRFv2_single_hemi_dif(CIy1,CIy2,LAI,lidf);  %the same for any
                
                ind=ind+1;
                BRFv2_hemi_dif_LUT(ind,:)=[sob_vsla_hemi_dif sof_vsla_hemi_dif kgd_dif];
            end
        end
    end
end

save BRFv2_hemi_dif_LUT.mat BRFv2_hemi_dif_LUT;



% %PART4
% id_LUT=zeros(6*6*70*8,1);
% 
% ind=0;
% for i_CIy1=1:6
%     CIy1=(i_CIy1-1)*0.1+0.5;
%     for i_CIy2=1:6
%         CIy2=(i_CIy2-1)*0.1+0.5;
%         for i_LAI=1:70
%             LAI=i_LAI*0.1;
%             for i_LIDFa=1:8
%                 LIDFa=i_LIDFa*10;
%                 lidf=campbell(LIDFa);
%                 lidf=lidf';
%                 
%                 id=i_hemi(CIy1,CIy2,LAI,lidf);  %the same for any
%                 
%                 ind=ind+1;
%                 id_LUT(ind,:)=id;
%             end
%         end
%     end
% end
% 
% save id_LUT.mat id_LUT;




