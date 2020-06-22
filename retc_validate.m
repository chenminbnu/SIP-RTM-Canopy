function LAIs = retc_validate
images = geotiffread('D:\科研\植被结构参数反演\data\retc\TM_20140616_clip.tif');
images = double(images)/10000;
landcover = geotiffread('D:\科研\植被结构参数反演\data\retc\landcover.tif');
load('RSO_LUT_Forest.mat');
load('RSO_LUT_Cropland.mat');

NIRs = images(:,:,5);
Reds = images(:,:,4);
LAIs = zeros(706,693);
for i = 1:706
%    tic
    for j = 1:693
        if(landcover(i,j)==1)
            LAIs(i,j)= RETC_Retrieval(Reds(i,j), NIRs(i,j),RSO_LUT_Forest);
        else
            LAIs(i,j)= RETC_Retrieval(Reds(i,j), NIRs(i,j),RSO_LUT_Cropland);
        end
    end
%    toc
end
end
