function [sob_vsla, sof_vsla,kg_dif]=A_BRFv2_single_dif(tto,CIo,CIy1,CIy2,LAI,lidf)

xx=[0.9602898565 -0.9602898565 0.7966664774 -0.7966664774 0.5255324099 -0.5255324099 0.1834346425 -0.1834346425];

ww=[0.1012285363 0.1012285363 0.2223810345 0.2223810345 0.3137066459 0.3137066459 0.3626837834 0.3626837834];

% * define limits of integration and the convertion factors for integration
% * over thetaL (note the tL suffix!)
upperlimit_tL = pi/2.0;
lowerlimit_tL = 0.0;
conv1_tL = (upperlimit_tL-lowerlimit_tL)/2.0;
conv2_tL = (upperlimit_tL+lowerlimit_tL)/2.0;

%   * define limits of integration and the convertion factors for integration
% * over phiL (note the pL suffix!)
upperlimit_pL = 2.0*pi;
lowerlimit_pL = 0.0;
conv1_pL = (upperlimit_pL-lowerlimit_pL)/2.0;
conv2_pL = (upperlimit_pL+lowerlimit_pL)/2.0;

sum_tL = zeros(1,1);
sum_tL_f = zeros(1,1);
sum_tL_g = zeros(1,1);

for i=1:8
    
    neword_tL = conv1_tL*xx(i) + conv2_tL;
    mu_tL     = cos(neword_tL);
    sin_tL    = sin(neword_tL);
    
    sum_pL = zeros(1,1);
    sum_pL_f = zeros(1,1);
    sum_pL_g = zeros(1,1);

    
    for j=1:8
        
        neword_pL  = conv1_pL*xx(j) + conv2_pL;
        
        tta=neword_tL*180/pi;              % observer zenith angle
        psia=neword_pL*180/pi;              % relative azimuth angle
        
        [Ga,Go,ka,K,sob,sof]    =   PROSAIL(tta,tto,psia,lidf);          
        
        [CIa]=CIxy(CIy1,CIy2,tta);
        
        [kca, kga]    =   sunshade(tta,tto,psia,Ga,Go,CIa,CIo,LAI);     
        
        sum_pL     = sum_pL + ww(j)*sob.*kca/K/pi;
        sum_pL_f     = sum_pL_f + ww(j)*sof.*kca/K/pi;
        sum_pL_g     = sum_pL_g + ww(j)*kga/pi;
    end
    
    sum_pL = sum_pL*conv1_pL;
          
    sum_tL = sum_tL + ww(i)* mu_tL*sin_tL*sum_pL;
    
    sum_pL_f = sum_pL_f*conv1_pL;
    
    sum_tL_f = sum_tL_f + ww(i)* mu_tL*sin_tL*sum_pL_f;
    
    sum_pL_g = sum_pL_g*conv1_pL;
    
    sum_tL_g = sum_tL_g + ww(i)* mu_tL*sin_tL*sum_pL_g;
    

end
sum_tL = sum_tL*conv1_tL;
sum_tL_f = sum_tL_f*conv1_tL;
sum_tL_g = sum_tL_g*conv1_tL;

sob_vsla=sum_tL;
sof_vsla=sum_tL_f;
kg_dif=sum_tL_g;
