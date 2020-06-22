%pp&cp
va=zeros(29*2,4);
va2=zeros(29*2,5);
for t=1:15
    va(t,1)=5*(15-t);
    va(t,2)=0+0;
end

for t=1:14
    va(t+15,1)=5*t;
    va(t+15,2)=0+180;
end
for t=1:15
    va(t+29,1)=5*(15-t);
    va(t+29,2)=0+270;
end

for t=1:14
    va(t+15+29,1)=5*t;
    va(t+15+29,2)=90;
end

for t=1:29*2
    
    tts		=	30.;		% solar zenith angle (?
    tto		=	va(t,1);		% observer zenith angle (?
    psi		=	va(t,2);         % azimuth (?
    
    if psi>180
        psi=psi-360;
    end
    psi=abs(psi);
    
    % % direct illumination
    % tts=0;
    % tto=0;
    % psi=0;
    LIDFa=30;
    rho=[0.0706 0.5756];
    tau=[0.0119 0.3379];
    LAI=3;
    
    w=rho+tau;
    
    rg=[0.3071 0.4102];
    
    CIy1=1;
    CIy2=1;
    
    [CIs]=CIxy(CIy1,CIy2,tts);
    [CIo]=CIxy(CIy1,CIy2,tto);
    
    deg2rad = pi/180;
    
    mius     = cos(tts*deg2rad);             %           cos solar       angle
    miuo     = cos(tto*deg2rad);             %           cos observation angle
    
    sins     = sin(tts*deg2rad);             %           cos solar       angle
    sino     = sin(tto*deg2rad);             %           cos observation angle
    
    psi         = abs(psi-360*round(psi/360));  %           (to ensure that volscatt is symmetric for psi=90 and psi=270)
    
    lidf=campbell(LIDFa);
    lidf=lidf';
    [Gs,Go,k,K,sob,sof]    =   PROSAIL(tts,tto,psi,lidf);
      
    i0=1-exp(-k*LAI*CIs);
    iv=1-exp(-K*LAI*CIo);
    
    t0=1-i0;
    tv=1-iv;
    
    [kc, kg]    =   sunshade(tts,tto,psi,Gs,Go,CIs,CIo,LAI);
    
    [sob_vsla, sof_vsla,kgd]=A_BRFv2_single_hemi(tts,CIs,CIy1,CIy2,LAI,lidf);  %the same for VZA    
    
    [sob_vsla_dif, sof_vsla_dif, kg_dif]=A_BRFv2_single_dif(tto,CIo,CIy1,CIy2,LAI,lidf);  %the same for VZA  
    
    [sob_vsla_hemi_dif, sof_vsla_hemi_dif,kgd_dif]=A_BRFv2_single_hemi_dif(CIy1,CIy2,LAI,lidf);  %the same for any  
    
    %direct illumination
    %so
%     SLA=i0/k;
    
%     rho1=kc/2/SLA;
    rho2=iv/2/LAI;
    
    id=i_hemi(CIy1,CIy2,LAI,lidf);  %the same for any
    td=1-id;
    
    p=1-id/LAI;  %its fixed
    
    %sd
%     rho_hemi1=kcd/2/SLA;
    rho_hemi2=id/2/LAI;
    
    
    %diffuse illumination
    %do
%     SLA_dif=sla_dif(CIy1,CIy2,LAI,lidf);  %be careful to use
    
%     rho_dif1=kc_dif/2/SLA_dif;
    rho_dif2=iv/2/LAI;
    
    %dd
%     rho_dif_hemi1=kcd_dif/2/SLA_dif;
    rho_dif_hemi2=id/2/LAI; %change from iv to id
    
    
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
    
    %sd
%     Rv=i0*w.*(rho_hemi1+p*w*rho_hemi2./(1-p*w));
    
    Rv2=sob_vsla*rho + sof_vsla*tau+i0*w.^2*p*rho_hemi2./(1-p*w); %test  wso.*kc/K    wso           = sob*rho + sof*tau;
      
    Rs=kgd*rg;
    
    Tup_hemi=td+id*w*rho_hemi2./(1-p*w);
    
    Rm=rg.*Tdn.*Tup_hemi./(1-rg.*Rdn)-t0*rg*td;
    
%     R=Rv+Rs+Rm;
    R2=Rv2+Rs+Rm;
    
    %absorption
    Av=i0*(1-w)./(1-p*w);
    Aup=id*(1-w)./(1-p*w);
    Am=rg.*Tdn.*Aup./(1-rg.*Rdn);
    A=Av+Am;
    
    
    %%results 2
    %diffuse illumination
    %do
%     BRF_difv=id*w.*(rho_dif1+p*w*rho_dif2./(1-p*w));
    
    BRF_difv2=sob_vsla_dif*rho + sof_vsla_dif*tau+id*w.^2*p*rho_dif2./(1-p*w);      
    
    BRF_difs=kg_dif*rg;
    
    Tdn_dif=td+id*w*rho_dif_hemi2./(1-p*w);
    Tup_difo=tv+id*w*rho_dif2./(1-p*w);
    Rdn_dif=id*w*rho_dif_hemi2./(1-p*w);
    
    BRF_difm=rg.*Tdn_dif.*Tup_difo./(1-rg.*Rdn_dif)-td*rg*tv;
    
%     BRF_dif=BRF_difv+BRF_difs+BRF_difm;
    
    BRF_dif2=BRF_difv2+BRF_difs+BRF_difm;
    
    %dd
%     R_difv=id*w.*(rho_dif_hemi1+p*w*rho_dif_hemi2./(1-p*w));
    
    R_difv2=sob_vsla_hemi_dif*rho + sof_vsla_hemi_dif*tau+id*w.^2*p*rho_dif_hemi2./(1-p*w);
    
    R_difs=kgd_dif*rg;
    
    Tup_dif_hemi=td+id*w*rho_dif_hemi2./(1-p*w);
    
    R_difm=rg.*Tdn_dif.*Tup_dif_hemi./(1-rg.*Rdn_dif)-td*rg*td;
    
%     R_dif=R_difv+R_difs+R_difm;
    
    R_dif2=R_difv2+R_difs+R_difm;
    
    %absorption
    A_difv=id*(1-w)./(1-p*w);
    Aup_dif=id*(1-w)./(1-p*w);
    A_difm=rg.*Tdn_dif.*Aup_dif./(1-rg.*Rdn_dif);
    A_dif=A_difv+A_difm;
    
    va2(t,1:2)=BRF;
    va2(t,3:4)=BRF2;
    
    va2(t,5)=kc;
end

% va(:,3:4)=va2(:,1:2);

va(:,3:4)=va2(:,3:4);

x=zeros(29,1);
for t=1:15
    x(t,1)=-5*(15-t);
end

for t=1:14
    x(t+15,1)=5*t;
end

subplot(2,2,1),plot(x,va(1:29,3),'r','Linewidth',2);
subplot(2,2,2),plot(x,va(29+1:29+29,3),'k','Linewidth',2);
subplot(2,2,3),plot(x,va(1:29,4),'r','Linewidth',2);
subplot(2,2,4),plot(x,va(29+1:29+29,4),'k','Linewidth',2);

