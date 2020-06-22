LAIs = retc_validate;
images = geotiffread('D:\科研\植被结构参数反演\data\retc\TM_20140616_clip.tif');
images = double(images)/10000;
landcover = geotiffread('D:\科研\植被结构参数反演\data\retc\landcover.tif');
load('RSO_LUT_Forest.mat');
load('RSO_LUT_Cropland.mat');

NIRs = images(:,:,5);
Reds = images(:,:,4);
NIR_simulation = zeros(706,693);
Red_simulation = zeros(706,693);
for i = 1:706
    tic
    for j = 1:693
        tmp = round(LAIs(i,j)*10);
        if tmp<=0
            tmp = 1;
        elseif tmp>70
            tmp = 70;
        end
        if(landcover(i,j)==1)
            RSO_LUT = RSO_LUT_Forest;
        else
           RSO_LUT = RSO_LUT_Cropland;
        end
        Red_simulation(i,j)= RSO_LUT(tmp,1);
        NIR_simulation(i,j)= RSO_LUT(tmp,2);
    end
    toc
end
% NIRs(landcover==1) = 0;
% Reds(landcover==1) = 0;
% NIR_simulation(landcover==1) = 0;
% Red_simulation(landcover==1) = 0;

%% plot
figure; 
subplot(2,3,1); imagesc(Reds,[0 0.06])
subplot(2,3,2); imagesc(Red_simulation,[0 0.06])
subplot(2,3,3); imagesc(Reds-Red_simulation,[-0.06 0.06])
subplot(2,3,4); imagesc(NIRs,[0 0.6])
subplot(2,3,5); imagesc(NIR_simulation,[0 0.6])
subplot(2,3,6); imagesc(NIRs-NIR_simulation,[-0.2 0.2])
figure; 
subplot(2,2,1);
hist(Reds(landcover==1)-Red_simulation(landcover==1),50)
subplot(2,2,2);
hist(NIRs(landcover==1)-NIR_simulation(landcover==1),50)
subplot(2,2,3);
hist(Reds(landcover==2)-Red_simulation(landcover==2),50)
subplot(2,2,4);
hist(NIRs(landcover==2)-NIR_simulation(landcover==2),50)