
clear;
rho=[0.0726 0.4409];
tau=[0.0726 0.4409];

rgs=[0.015000000  0.048185132
    0.075000003  0.114972815
    0.119999997  0.165063575
    0.165000007  0.215154335
    0.209999993  0.265245110
    0.254999995  0.315335870
    0.300000012  0.365426630
    0.344999999  0.415517390
    0.389999986  0.465608150
    0.449999988  0.532395840];
w=rho+tau;

load BRFv2_dif_LUT.mat;
load BRFv2_hemi_LUT.mat;
load BRFv2_hemi_dif_LUT.mat;
load id_LUT.mat;

%%RSO
RSO_LUT=zeros(7*7*7*6*6*70*8*10,6);
ind=0;
%%parpool(24)
%par
for i_tts=1:7
    tts=(i_tts-1)*10+5;
    for i_tto=1:7
        tto=(i_tto-1)*10;
        for i_psi=1:7
            psi=(i_psi-1)*30;
            tic
            for i_CIy1=1:6
                for i_CIy2=1:6
                    CIy1=(i_CIy1-1)*0.1+0.5;
                    CIy2=(i_CIy2-1)*0.1+0.5;
                    [CIs]=CIxy(CIy1,CIy2,tts);
                    [CIo]=CIxy(CIy1,CIy2,tto);
                    for i_LAI=1:70
                        LAI=i_LAI*0.1;
                        for i_LIDFa=1:8
                            LIDFa=i_LIDFa*10;
                            lidf=campbell(LIDFa);
                            lidf=lidf';
                            for i_rsoil = 1:10
                                try
                                    rg = rgs(i_rsoil,:);
                                    %                             psi         = abs(psi-360*round(psi/360));  %           (to ensure that volscatt is symmetric for psi=90 and psi=270)
                                    [Gs,Go,k,K,sob,sof]    =   PROSAIL(tts,tto,psi,lidf);
                                    i0=1-exp(-k*LAI*CIs);
                                    iv=1-exp(-K*LAI*CIo);
                                    t0=1-i0;
                                    tv=1-iv;
                                    [kc, kg]    =   sunshade(tts,tto,psi,Gs,Go,CIs,CIo,LAI);
                                    %direct illumination
                                    %so
                                    %     SLA=i0/k;
                                    %     rho1=kc/2/SLA;
                                    rho2=iv/2/LAI;
                                    %     id=i_hemi(CIy1,CIy2,LAI,lidf);  %the same for any
                                    id=id_LUT((i_CIy1-1)*6*70*8+(i_CIy2-1)*70*8+(i_LAI-1)*8+i_LIDFa);  %the same for any
                                    td=1-id;
                                    p=1-id/LAI;  %its fixed
                                    
                                    %sd
                                    %     rho_hemi1=kcd/2/SLA;
                                    rho_hemi2=id/2/LAI;
                                    
                                    %%results
                                    %direct illumination
                                    %so
                                    wso           = sob*rho + sof*tau;            % [nwl]     w,      p309{1} bidirectional scattering coefficent (directional-directional)
                                    %     BRFv=i0*w.*(rho1+p*w*rho2./(1-p*w));
                                    BRFv2=wso.*kc/K+i0*w.^2*p*rho2./(1-p*w);   %this is more accurate
                                    BRFs=kg*rg;
                                    Tdn=t0+i0*w*rho_hemi2./(1-p*w);
                                    Tup_o=tv+id*w*rho2./(1-p*w);
                                    Rdn=id*w*rho_hemi2./(1-p*w);
                                    BRFm=rg.*Tdn.*Tup_o./(1-rg.*Rdn)-t0*rg*tv;
                                    %     BRF=BRFv+BRFs+BRFm;
                                    BRF2=BRFv2+BRFs+BRFm;
                                    
                                    %absorption
                                    Av=i0*(1-w)./(1-p*w);
                                    Aup=id*(1-w)./(1-p*w);
                                    Am=rg.*Tdn.*Aup./(1-rg.*Rdn);
                                    A=Av+Am;
                                    
                                    ind=ind+1;
                                    RSO_LUT(ind,1:5)=[BRF2 A i0];
                                catch exception
                                    disp(exception);
                                end
                            end
                        end
                    end
                end
                
            end
            toc
        end
    end
end

RSO_LUT(:,6)=(RSO_LUT(:,2)-RSO_LUT(:,1))./(RSO_LUT(:,2)+RSO_LUT(:,1)).*RSO_LUT(:,2);
save('RSO_LUT_update.mat','RSO_LUT','-v7.3');


