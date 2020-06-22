function Retrieval_LAI(index)
%% Read EPIC MAIAC surface reflectance, AOD, and Kernel Parameters

%% constant
cd('/pic/projects/sif/dalei/code/SIP');
load('RSO_LUT_update.mat');
%clc; clear all;
%% input
hi = index;
Tm={'Tile00','Tile01','Tile02','Tile03','Tile10','Tile11','Tile12','Tile13'};
tm={'tile00','tile01','tile02','tile03','tile10','tile11','tile12','tile13'};

%for hi=1:8 % tile
disp(hi)
%% Output
filepath=['/pic/projects/sif/dalei/EPIC-data/Extract_Reflectance/',Tm{hi},'/'];
fileOutpath=['/pic/projects/sif/dalei/EPIC-data/Estimate LAI/update/',Tm{hi},'/'];

for t = [185:5:365]
    
    disp(t)
    TimeStart = datetime(2016,1,1,0,0,0);
    TimeStart = TimeStart + t - 5;
    NIRs = [];
    Reds = [];
    SZAs = [];
    VZAs = [];
    RAAs = [];
    index = 0;
    for day= 1:5
    disp(day);
        IndexTime = TimeStart + day-1;
        [Indexyear,IndexMonth, Indexday] = ymd(IndexTime);
        filestr = dir([filepath '*' num2str(Indexyear) num2str(IndexMonth,'%02d')  num2str(Indexday,'%02d')   '*.mat']);
  % disp([filepath '*' num2str(Indexyear) num2str(IndexMonth,'%02d')  num2str(Indexday,'%02d')   '*.mat']);
         [filenum,nouse]=size(filestr );
        if filenum>0
            for fi=1:filenum
                load([filepath filestr(fi).name]);
                filename_tmp = filestr(fi).name;
                
                
                filter_SZA = isnan(SZA1);
                SZA_tmp = SZA1*180/pi;
                SZA_tmp(SZA1>89)=nan;
                SZA_tmp(SZA1<0)=nan;
                SZA_tmp(filter_SZA) = nan;

                filter_VZA = isnan(VZA1);
                VZA_tmp = VZA1*180/pi;
                VZA_tmp(VZA1>89)=nan;
                VZA_tmp(VZA1<0)=nan;
                VZA_tmp(filter_VZA) = nan;

                RAA_tmp = abs(SAA1-VAA1)*180/pi;
                %RAA_tmp(RAA_tmp>180) = 360-RAA_tmp(RAA_tmp>180);
                filter_RAA = isnan(SAA1) |isnan(VAA1) ;
                RAA_tmp(RAA_tmp>360)=nan;
                RAA_tmp(RAA_tmp<0)=nan;
                RAA_tmp(filter_RAA) = nan;
                NIR_tmp = NIR;
                Red_tmp = Red;
                %% calculate SW and PAR albedo
                filter_Q =  ~((Cloud_mask==1)  & (LW_mask>=1) & (LW_mask<=5) & (Status_QA==0));
                
                Red_tmp(filter_Q) = nan;
                NIR_tmp(filter_Q) = nan;
                RAA_tmp(filter_Q) = nan;
                VZA_tmp(filter_Q) = nan;
                SZA_tmp(filter_Q) = nan;
                index = index + 1;
                Reds(:,:,index) = Red_tmp;
                NIRs(:,:,index) = NIR_tmp;
                RAAs(:,:,index) = RAA_tmp;
                VZAs(:,:,index) = VZA_tmp;
                SZAs(:,:,index) = SZA_tmp;
            end
        end
        
    end
disp('start retrieving');
tic
    %% retrieval pixel by pixel
    rows = size(Reds,1);
    cols = size(NIRs,2);
    Results = zeros(rows, cols, 6);
    for row = 1:rows
        for col = 1:cols
            
              Red_tmp= squeeze(Reds(row,col,:));
                NIR_tmp=squeeze(NIRs(row,col,:));
                SZA_tmp = squeeze(SZAs(row,col,:));
                VZA_tmp = squeeze(VZAs(row,col,:));
                RAA_tmp = squeeze(RAAs(row,col,:));
            filter = (~isnan( Red_tmp)) & (~isnan( NIR_tmp))& (~isnan( SZA_tmp))& (~isnan( VZA_tmp))& (~isnan( RAA_tmp) & (SZA_tmp<=55)& (VZA_tmp<=55));
         if(sum(filter)>0)
              Red_tmp= Red_tmp(filter);
                 NIR_tmp=NIR_tmp(filter);
                SZA_tmp=SZA_tmp(filter);
                VZA_tmp= VZA_tmp(filter);
                RAA_tmp= RAA_tmp(filter);

            retrieval_value = A2_Retrieval_LUT_Step3(Red_tmp,NIR_tmp,SZA_tmp, VZA_tmp , RAA_tmp,RSO_LUT);
            Results(row,col,:) = [sum(filter);retrieval_value];
         end        
end
    end

toc
    save([fileOutpath 'red_nir_alter_2016-' num2str(hi,'%02d') '-d' num2str(t,'%03d') '.mat'],'Results');
end
disp('finished!');

end

