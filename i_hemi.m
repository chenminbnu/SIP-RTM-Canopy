function [sum_tL]=i_hemi(CIy1,CIy2,LAI,lidf)

xx=[0.9602898565 -0.9602898565 0.7966664774 -0.7966664774 0.5255324099 -0.5255324099 0.1834346425 -0.1834346425];

ww=[0.1012285363 0.1012285363 0.2223810345 0.2223810345 0.3137066459 0.3137066459 0.3626837834 0.3626837834];

% * define limits of integration and the convertion factors for integration
% * over thetaL (note the tL suffix!)
upperlimit_tL = pi/2.0;
lowerlimit_tL = 0.0;
conv1_tL = (upperlimit_tL-lowerlimit_tL)/2.0;
conv2_tL = (upperlimit_tL+lowerlimit_tL)/2.0;

% %   * define limits of integration and the convertion factors for integration
% % * over phiL (note the pL suffix!)
% upperlimit_pL = 2.0*pi;
% lowerlimit_pL = 0.0;
% conv1_pL = (upperlimit_pL-lowerlimit_pL)/2.0;
% conv2_pL = (upperlimit_pL+lowerlimit_pL)/2.0;

sum_tL = zeros(1,1);

for i=1:8
    
    neword_tL = conv1_tL*xx(i) + conv2_tL;
    mu_tL     = cos(neword_tL);
    sin_tL    = sin(neword_tL);
    
%     sum_pL = zeros(1,1);

    
%     for j=1:8
        
%         neword_pL  = conv1_pL*xx(j) + conv2_pL;
        
        tta=neword_tL*180/pi;              % observer zenith angle
%         psia=neword_pL*180/pi;              % relative azimuth angle
        
        [Ga,ka]    =   PROSAIL2(tta,lidf); 
        
        [CIa]=CIxy(CIy1,CIy2,tta);
        
        ia=1-exp(-ka*LAI*CIa);
               
%         sum_pL     = sum_pL + ww(j)*ia/pi;
%     end
    
%     sum_pL = sum_pL*conv1_pL;
          
%     sum_tL = sum_tL + ww(i)* mu_tL*sin_tL*sum_pL;

    sum_tL = sum_tL + ww(i)* mu_tL*sin_tL*ia*2;

end
sum_tL = sum_tL*conv1_tL;
