LAIs = retc_validate;
images = geotiffread('D:\科研\植被结构参数反演\data\retc\TM_20140616_clip.tif');
images = double(images)/10000;
landcover = geotiffread('D:\科研\植被结构参数反演\data\retc\landcover.tif');
 load('RSO_LUT_Forest.mat');
 load('RSO_LUT_Cropland.mat');
NIRs = images(:,:,5);
Reds = images(:,:,4);
NIRs_coarse = zeros(22,21);
REDs_coarse = zeros(22,21);
NIRs_coarse2 = zeros(22,21);
REDs_coarse2 = zeros(22,21);

NIRs_model1 =  zeros(22,21);
NIRs_model2 =  zeros(22,21);
REDs_model1 =  zeros(22,21);
REDs_model2 =  zeros(22,21);
Landcover_coarse = zeros(22,21);
LAIs_coarse1 = zeros(22,21);
LAIs_coarse2_f = zeros(22,21);
LAIs_coarse2_c = zeros(22,21);
for i = 1:44
    for j = 1:42
        NIRs_coarse(i,j) = mean(mean(NIRs((i-1)*16+1:i*16,(j-1)*16+1:j*16)));
        REDs_coarse(i,j) = mean(mean(Reds((i-1)*16+1:i*16,(j-1)*16+1:j*16)));
           NIRs_coarse2(i,j) = mean(mean(NIR_simulation((i-1)*16+1:i*16,(j-1)*16+1:j*16)));
        REDs_coarse2(i,j) = mean(mean(Red_simulation((i-1)*16+1:i*16,(j-1)*16+1:j*16)));
        Landcover_coarse(i,j) =  sum(sum(landcover((i-1)*16+1:i*16,(j-1)*16+1:j*16)==1))/(16*16);
        tmp = LAIs((i-1)*16+1:i*16,(j-1)*16+1:j*16);
        tmp_class = landcover((i-1)*16+1:i*16,(j-1)*16+1:j*16);
        LAIs_coarse1(i,j) = mean(mean(tmp));
        LAIs_coarse2_f(i,j) = mean(tmp(tmp_class==1));
        LAIs_coarse2_c(i,j) = mean(tmp(tmp_class==2));
      if(isnan( LAIs_coarse2_f(i,j)))
          LAIs_coarse2_f(i,j) = 0;
      end
      if(isnan( LAIs_coarse2_c(i,j)))
          LAIs_coarse2_c(i,j) = 0;
      end
        
        %% forward simulation 
      
       % simulation 1
       tmp = round( LAIs_coarse1(i,j)*10);
        if tmp<=0
            tmp = 1;
        elseif tmp>70
            tmp = 70;
        end
        if( Landcover_coarse(i,j)>=0.5)
            RSO_LUT = RSO_LUT_Forest;
        else
           RSO_LUT = RSO_LUT_Cropland;
        end
        NIRs_model1(i,j) =  RSO_LUT(tmp,2);
        REDs_model1(i,j) =  RSO_LUT(tmp,1);
          %% simulation 2
         tmp_1 = round( LAIs_coarse2_f(i,j)*10);
           tmp_2 = round( LAIs_coarse2_c(i,j)*10);
        if tmp_1<=0
            tmp_1 = 1;
        elseif tmp_1>70
            tmp_1 = 70;
        end
         if tmp_2<=0
            tmp_2 = 1;
        elseif tmp_2>70
            tmp_2= 70;
         end
          tmp_1-tmp_2
            RSO_LUT = RSO_LUT_Forest;
        tmp1_red= RSO_LUT(tmp_1,1);
        tmp1_nir= RSO_LUT(tmp_1,2);
         RSO_LUT = RSO_LUT_Cropland;
        tmp2_red= RSO_LUT(tmp_2,1);
        tmp2_nir= RSO_LUT(tmp_2,2);
        
        NIRs_model2(i,j)= tmp1_nir*Landcover_coarse(i,j) + tmp2_nir*(1-Landcover_coarse(i,j)) ;
        REDs_model2(i,j)= tmp1_red*Landcover_coarse(i,j) + tmp2_red*(1-Landcover_coarse(i,j));
    
        
    end
end
%% plot
figure;
subplot(2,3,1); imagesc(REDs_coarse,[0 0.06])
subplot(2,3,2); imagesc(REDs_model1,[0 0.06])
subplot(2,3,3); imagesc(REDs_model2,[0 0.06])
subplot(2,3,4); imagesc(NIRs_coarse,[0 0.6])
subplot(2,3,5); imagesc(NIRs_model1,[0 0.6])
subplot(2,3,6); imagesc(NIRs_model2,[0 0.6])
figure;
subplot(2,2,1); imagesc(REDs_coarse-REDs_model1,[-0.06 0.06])
subplot(2,2,2); imagesc(REDs_coarse-REDs_model2,[-0.06 0.06])
subplot(2,2,3); imagesc(NIRs_coarse-NIRs_model1,[-0.1 0.1])
subplot(2,2,4); imagesc(NIRs_coarse-NIRs_model2,[-0.1 0.1])

figure;
subplot(2,2,1); hist(REDs_coarse(:)-REDs_model1(:),50)
subplot(2,2,2); hist(REDs_coarse(:)-REDs_model2(:),50)
subplot(2,2,3); hist(NIRs_coarse(:)-NIRs_model1(:),50)
subplot(2,2,4); hist(NIRs_coarse(:)-NIRs_model2(:),50)

figure
subplot(1,2,1)
hold on
scatter(Landcover_coarse(Landcover_coarse>0),NIRs_coarse(Landcover_coarse>0)-NIRs_model1(Landcover_coarse>0),'r')
scatter(Landcover_coarse(Landcover_coarse>0),NIRs_coarse(Landcover_coarse>0)-NIRs_model2(Landcover_coarse>0),'b')
legend('ordinary','RTEC')
xlabel('percent of forest')
ylabel('error')

title('NIR')
subplot(1,2,2)
hold on
scatter(Landcover_coarse(:),REDs_coarse(:)-REDs_model1(:),'r')
scatter(Landcover_coarse(:),REDs_coarse(:)-REDs_model2(:),'b')
legend('ordinary','RTEC')
title('Red')
xlabel('percent of forest')
ylabel('error')

%% figure
figure
subplot(1,2,1)
hold on
scatter(Landcover_coarse(Landcover_coarse>0),NIRs_coarse2(Landcover_coarse>0)-NIRs_model1(Landcover_coarse>0),'r')
scatter(Landcover_coarse(Landcover_coarse>0),NIRs_coarse2(Landcover_coarse>0)-NIRs_model2(Landcover_coarse>0),'b')
legend('ordinary','RTEC')
xlabel('percent of forest')
ylabel('error')

title('NIR')
subplot(1,2,2)
hold on
scatter(Landcover_coarse(Landcover_coarse>0),REDs_coarse2(Landcover_coarse>0)-REDs_model1(Landcover_coarse>0),'r')
scatter(Landcover_coarse(Landcover_coarse>0),REDs_coarse2(Landcover_coarse>0)-REDs_model2(Landcover_coarse>0),'b')
legend('ordinary','RTEC')
title('Red')
xlabel('percent of forest')
ylabel('error')
%% figure
figure
subplot(1,2,1)
hold on
scatter(NIRs_coarse(Landcover_coarse>0),NIRs_model1(Landcover_coarse>0),'r')
refline
scatter(NIRs_coarse(Landcover_coarse>0),NIRs_model2(Landcover_coarse>0),'b')
refline
axis([0.2 0.5 0.2 0.5])
legend('ordinary','RTEC')
xlabel('percent of forest')
ylabel('error')
plot([0.2,0.5],[0.2 0.5])

title('NIR')
subplot(1,2,2)
hold on
scatter(REDs_coarse(Landcover_coarse>0),REDs_model1(Landcover_coarse>0),'r')
scatter(REDs_coarse(Landcover_coarse>0),REDs_model2(Landcover_coarse>0),'b')
axis([0 0.1 0 0.1])

legend('ordinary','RTEC')
title('Red')
xlabel('percent of forest')
ylabel('error')
plot([0,0.1],[0 0.1])

%% figure
t = Landcover_coarse(Landcover_coarse>=0.1 & Landcover_coarse<=0.9);
t(t<0.5) = 1-t(t<0.5);
figure
colormap jet
subplot(2,2,2)
hold on
plot([0,0.1],[0 0.1],'k','Linewidth',1.5)

scatter(REDs_coarse2(Landcover_coarse>=0.1 & Landcover_coarse<=0.9),REDs_model1(Landcover_coarse>=0.1 & Landcover_coarse<=0.9),15,t,'filled')
t1 = REDs_coarse2(Landcover_coarse>=0.1 & Landcover_coarse<=0.9);

t2 = REDs_model1(Landcover_coarse>=0.1 & Landcover_coarse<=0.9);
R2 = corrcoef(REDs_coarse2(Landcover_coarse>=0.1 & Landcover_coarse<=0.9),REDs_model1(Landcover_coarse>=0.1 & Landcover_coarse<=0.9));
R2 = R2(1,2);
R2  = R2^2;
RMSE = sqrt(mean((t1-t2).^2));
text(0.005,0.09,['R^2=' num2str(R2,'%.2f')])
text(0.005,0.078,['RMSE=' num2str(RMSE,'%.3f')])
axis([0 0.1 0 0.1])
set(gca,'linewidth',1.5)
%legend('ordinary','RTEC')
title('Red')
xlabel('Reference Reflectance')
ylabel('SRTE Reflectance')
box on
colorbar
caxis([0.5 1.0])
subplot(2,2,1)
hold on
plot([0,0.1],[0 0.1],'k','Linewidth',1.5)

scatter(REDs_coarse2(Landcover_coarse>=0.1 & Landcover_coarse<=0.9),REDs_model2(Landcover_coarse>=0.1 & Landcover_coarse<=0.9),15,t,'filled')
axis([0 0.1 0 0.1])
t1 = REDs_coarse2(Landcover_coarse>=0.1 & Landcover_coarse<=0.9);
t2 = REDs_model2(Landcover_coarse>=0.1 & Landcover_coarse<=0.9);
R2 = corrcoef(REDs_coarse2(Landcover_coarse>=0.1 & Landcover_coarse<=0.9),REDs_model2(Landcover_coarse>=0.1 & Landcover_coarse<=0.9));
R2 = R2(1,2);
R2  = R2^2;
RMSE = sqrt(mean((t1-t2).^2));
text(0.005,0.09,['R^2=' num2str(R2,'%.2f')])
text(0.005,0.078,['RMSE=' num2str(RMSE,'%.3f')])
set(gca,'linewidth',1.5)

%legend('ordinary','RTEC')
title('Red')
xlabel('Reference Reflectance')
ylabel('RTEC Reflectance')
box on
colorbar
caxis([0.5 1.0])
subplot(2,2,4)
hold on
plot([0.3,0.5],[0.3 0.5],'k','Linewidth',1.5)
scatter(NIRs_coarse2(Landcover_coarse>=0.1 & Landcover_coarse<=0.9),NIRs_model1(Landcover_coarse>=0.1 & Landcover_coarse<=0.9),15,t,'filled')
t1 = NIRs_coarse2(Landcover_coarse>=0.1 & Landcover_coarse<=0.9);
t2 = NIRs_model1(Landcover_coarse>=0.1 & Landcover_coarse<=0.9);
R2 = corrcoef(NIRs_coarse2(Landcover_coarse>=0.1 & Landcover_coarse<=0.9),NIRs_model1(Landcover_coarse>=0.1 & Landcover_coarse<=0.9));
R2 = R2(1,2);
R2  = R2^2;
RMSE = sqrt(mean((t1-t2).^2));
text(0.41,0.345,['R^2=' num2str(R2,'%.2f')])
text(0.41,0.32,['RMSE=' num2str(RMSE,'%.3f')])
set(gca,'linewidth',1.5)

axis([0.3 0.5 0.3 0.5])
%legend('ordinary','RTEC')
xlabel('Reference Reflectance')
ylabel('SRTE Reflectance')
box on
%title('NIR')
colorbar
caxis([0.5 1.0])
subplot(2,2,3)
% refline
hold on
plot([0.3,0.5],[0.3 0.5],'k','Linewidth',1.5)

scatter(NIRs_coarse2(Landcover_coarse>=0.1 & Landcover_coarse<=0.9),NIRs_model2(Landcover_coarse>=0.1 & Landcover_coarse<=0.9),15,t,'filled')
t1 = NIRs_coarse2(Landcover_coarse>=0.1 & Landcover_coarse<=0.9);
t2 = NIRs_model2(Landcover_coarse>=0.1 & Landcover_coarse<=0.9);
R2 = corrcoef(NIRs_coarse2(Landcover_coarse>=0.1 & Landcover_coarse<=0.9),NIRs_model2(Landcover_coarse>=0.1 & Landcover_coarse<=0.9));
R2 = R2(1,2);
R2  = R2^2;
RMSE = sqrt(mean((t1-t2).^2));
text(0.41,0.345,['R^2=' num2str(R2,'%.2f')])
text(0.41,0.32,['RMSE=' num2str(RMSE,'%.3f')])

% refline
axis([0.3 0.5 0.3 0.5])
%legend('ordinary','RTEC')
xlabel('Reference Reflectance')
ylabel('RTEC Reflectance')
box on
title('NIR')
set(gca,'linewidth',1.5)
colorbar
caxis([0.5 1.0])