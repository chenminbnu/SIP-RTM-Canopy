function [kc, kg]    =   sunshade(tts,tto,psi,Gs,Go,CIs,CIo,LAI)

deg2rad = pi/180;

% mius     = cos(tts*deg2rad);             %           cos solar       angle 
% miuo     = cos(tto*deg2rad);             %           cos observation angle 
% 
% sins     = sin(tts*deg2rad);             %           cos solar       angle 
% sino     = sin(tto*deg2rad);             %           cos observation angle 
% 
% alp=acos(mius*miuo+sins*sino*cos(psi*deg2rad)); % 视线和光线间的夹角
% 
% siga=1/mius^2+1/miuo^2-2*cos(alp)/mius/miuo; % 为防止sig==0
% if siga<0.0000001
%     sig=0.0000001;
% else
%     sig=siga^0.5;
% end
% 
% d=0.05; % 与Kuusk文章中的dL含义相同?????????与叶片形状、叶倾角分布有关的可调参数
% H=1;
% W=d/H/sig*(1-exp(-H*sig/d));

% kc1=1-exp(-Gs*CIs/mius*LAI)-exp(-Go*CIo/miuo*LAI)+exp((-Gs*CIs/mius-Go*CIo/miuo+(Gs*Go*CIs*CIo/mius/miuo)^0.5*W)*LAI); % by zeng yelu计算阳叶面积比例
% kc2=1-exp((-(Gs*Go*CIs*CIo/mius/miuo)^0.5*W)*LAI); % by fan wenjie计算阳叶面积比例
% kc=(kc1+kc2)/2;  %取平均
% 
% kg=exp((-Gs*CIs/mius-Go*CIo/miuo+(Gs*Go*CIs*CIo/mius/miuo)^0.5*W)*LAI);

cos_tts     = cos(tts*deg2rad);             %           cos solar       angle   
tan_tto     = tan(tto*deg2rad);             %           tan observation angle

cos_tto     = cos(tto*deg2rad);             %           cos observation angle   
tan_tts     = tan(tts*deg2rad);             %           tan observation angle

psi         = abs(psi-360*round(psi/360));  %           (to ensure that volscatt is symmetric for psi=90 and psi=270)

if tts==tto&&psi==0
    kc=1-exp(-Gs*CIs/cos_tts*LAI);
    kg=exp(-Gs*CIs/cos_tts*LAI);
else
    
nl              = 20;
x        = (-1/nl : -1/nl : -1)';         % a column vector
xl       = [0; x];                 % add top level
dx      = 1/nl;
iLAI    = LAI/nl;               % [1]          LAI of elementary layer
d=0.05; % 与Kuusk文章中的dL含义相同?????????与叶片形状、叶倾角分布有关的可调参数
H=1;
q=d/H;

dso         = sqrt(tan_tts.^2 + tan_tto.^2 - 2*tan_tts.*tan_tto.*cos(psi*deg2rad));

k=Gs/cos_tts;
K=Go/cos_tto;

Ps          =   exp(k*xl*CIs*LAI);                                              % [nl+1]  p154{1} probability of viewing a leaf in solar dir
Po          =   exp(K*xl*CIo*LAI);                                              % [nl+1]  p154{1} probability of viewing a leaf in observation dir

Ps(1:nl)    =   Ps(1:nl) *(1-exp(-k*CIs*LAI*dx))/(k*CIs*LAI*dx);                                      % Correct Ps/Po for finite dx
Po(1:nl)    =   Po(1:nl) *(1-exp(-K*CIo*LAI*dx))/(K*CIo*LAI*dx);  % Correct Ps/Po for finite dx

Pso         =   zeros(size(Po));
for j=1:length(xl)
    Pso(j,:)=   quad(@(y)Psofunction(K,k,CIs,CIo,LAI,q,dso,y),xl(j)-dx,xl(j))/dx; %#ok<FREMO>
end

Pso(Pso>Po)= min([Po(Pso>Po),Ps(Pso>Po)],[],2);    %takes care of rounding error
Pso(Pso>Ps)= min([Po(Pso>Ps),Ps(Pso>Ps)],[],2);    %takes care of rounding error

kc=iLAI*CIo*K*sum(Pso(1:nl));  %visable sunlit leaf
kg=Pso(nl+1);  %visable sunlit soil
end


function pso    =   Psofunction(K,k,CIs,CIo,LAI,q,dso,xl)
if dso~=0
    alf         =   (dso/q) *2/(k+K);
    pso         =   exp((K*CIo+k*CIs)*LAI*xl + sqrt(K*CIo*k*CIs)*LAI/(alf  )*(1-exp(xl*(alf  ))));% [nl+1]  factor for correlation of Ps and Po
else
    pso         =   exp((K*CIo+k*CIs)*LAI*xl - sqrt(K*CIo*k*CIs)*LAI*xl);% [nl+1]  factor for correlation of Ps and Po
end
