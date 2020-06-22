
clear;
rho=[0.0670 0.4713];
tau=[0.0398 0.4536];

rg=[0.1843  0.2842];
% rho=[0.0581 0.4839];
% tau=[0.0390 0.4809];
rg=[0.1843 0.2842];
w=rho+tau;

% load BRFv2_dif_LUT.mat;
% load BRFv2_hemi_LUT.mat;
% load BRFv2_hemi_dif_LUT.mat;
load id_LUT.mat;

%%RSO
RSO_LUT_Forest=zeros(70,6);
ind=0;
tts=30.0839;
tto=0;
psi=147.26471553;
i_CIy1=3;
i_CIy2=3;
CIy1=(i_CIy1-1)*0.1+0.5;
CIy2=(i_CIy2-1)*0.1+0.5;
[CIs]=CIxy(CIy1,CIy2,tts);
[CIo]=CIxy(CIy1,CIy2,tto);

i_LIDFa=6;
LIDFa=57.3;
lidf=campbell(LIDFa);
lidf=lidf';
for i_LAI=1:70
    LAI=i_LAI*0.1;
    try
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
        RSO_LUT_Forest(ind,1:5)=[BRF2 A i0];
    catch exception
        disp(exception);
    end
end

RSO_LUT_Forest(:,6)=(RSO_LUT_Forest(:,2)-RSO_LUT_Forest(:,1))./(RSO_LUT_Forest(:,2)+RSO_LUT_Forest(:,1)).*RSO_LUT_Forest(:,2);
save('RSO_LUT_Forest.mat','RSO_LUT_Forest','-v7.3');


