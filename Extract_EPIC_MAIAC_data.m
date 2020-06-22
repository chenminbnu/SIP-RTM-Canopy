%% Read EPIC MAIAC surface reflectance, AOD, and Kernel Parameters


%clc; clear all;
%% input

filename1='DSCOVR_EPIC_L2_MAIAC_01_2018';
filename2='*_02.h5';

% Loc=zeros(6,4);
% Loc(1,:)=[1,111,223,334];
% Loc(2,:)=[335,446,670,781];
% Loc(3,:)=[557,668,111,222];
% Loc(4,:)=[890,1000,223,334];
% Loc(5,:)=[668,779,1,112];
% Loc(6,:)=[557,668,780,891];
% Aname={'Amazon','North','Mississippi','Africa','South China','Heihe'};

Tm={'Tile00','Tile01','Tile02','Tile03','Tile10','Tile11','Tile12','Tile13'};
tm={'tile00','tile01','tile02','tile03','tile10','tile11','tile12','tile13'};
Angle_factor=1;
BRDF_factor=0.0001;
rd=acos(-1)/180.0;
for hi=1:8 % tile
%% Output
fileoutpath=['/pic/projects/sif/dalei/EPIC-data/Extract_Reflectance/2018/',Tm{hi},'/'];

filepath='/pic/projects/sif/dalei/EPIC-data/MAIAC/2018/';

    filestr=dir([filepath,filename1,filename2]);
    [filenum,nouse]=size(filestr);
    if filenum>0
   
    for fi=1:filenum
  disp(fi);
       
        %em_info=h5info([filepath,filestr(fi).name]);   % check metadata
        %for attribute name
        check_avail=h5readatt([filepath,filestr(fi).name],'/',[tm{hi},'_availability']);
       filename_tmp = filestr(fi).name;
        mi = str2num(filename_tmp(29:30));
        dayi = str2num(filename_tmp(31:32));
        if check_avail==1
           Dataset_name=['/',tm{hi},'/Sur_refl_680nm'];
           Redtp=h5read([filepath,filestr(fi).name],Dataset_name); 
           Red=double(Redtp')*BRDF_factor;
           Dataset_name=['/',tm{hi},'/Sur_refl_780nm'];
           NIRtp=h5read([filepath,filestr(fi).name],Dataset_name); 
           NIR=double(NIRtp')*BRDF_factor;
                     
           Dataset_name=['/',tm{hi},'/Sur_refl_vzen'];
           VZAtp=h5read([filepath,filestr(fi).name],Dataset_name); 
           VZA1=VZAtp'*Angle_factor*rd;  
           
           Dataset_name=['/',tm{hi},'/Sur_refl_szen'];
           SZAtp=h5read([filepath,filestr(fi).name],Dataset_name); 
           SZA1=SZAtp'*Angle_factor*rd;
           
           Dataset_name=['/',tm{hi},'/Sur_refl_vaz'];
           VAAtp=h5read([filepath,filestr(fi).name],Dataset_name); 
           VAA1=VAAtp'*Angle_factor*rd;  
           
           Dataset_name=['/',tm{hi},'/Sur_refl_saz'];
           SAAtp=h5read([filepath,filestr(fi).name],Dataset_name); 
           SAA1=SAAtp'*Angle_factor*rd;
           
           Dataset_name=['/',tm{hi},'/Aer_optical_depth'];
           AODtp=h5read([filepath,filestr(fi).name],Dataset_name); 
           AOD=double(AODtp')*BRDF_factor;

           Dataset_name=['/',tm{hi},'/Cloud_mask'];
           Cloud_mask=h5read([filepath,filestr(fi).name],Dataset_name); 
           Cloud_mask  = Cloud_mask';

           Dataset_name=['/',tm{hi},'/Kgeo_443nm'];
           Kgeo_443nm=h5read([filepath,filestr(fi).name],Dataset_name); 
           Kgeo_443nm=double(Kgeo_443nm')*BRDF_factor;
           Dataset_name=['/',tm{hi},'/Kiso_443nm'];
           Kiso_443nm=h5read([filepath,filestr(fi).name],Dataset_name); 
           Kiso_443nm=double(Kiso_443nm')*BRDF_factor;
           Dataset_name=['/',tm{hi},'/Kvol_443nm'];
           Kvol_443nm=h5read([filepath,filestr(fi).name],Dataset_name); 
           Kvol_443nm=double(Kvol_443nm')*BRDF_factor;

           Dataset_name=['/',tm{hi},'/Kgeo_551nm'];
           Kgeo_551nm=h5read([filepath,filestr(fi).name],Dataset_name); 
           Kgeo_551nm=double(Kgeo_551nm')*BRDF_factor;
           Dataset_name=['/',tm{hi},'/Kiso_551nm'];
           Kiso_551nm=h5read([filepath,filestr(fi).name],Dataset_name); 
           Kiso_551nm=double(Kiso_551nm')*BRDF_factor;
           Dataset_name=['/',tm{hi},'/Kvol_551nm'];
           Kvol_551nm=h5read([filepath,filestr(fi).name],Dataset_name); 
           Kvol_551nm=double(Kvol_551nm')*BRDF_factor;

           Dataset_name=['/',tm{hi},'/Kgeo_680nm'];
           Kgeo_680nm=h5read([filepath,filestr(fi).name],Dataset_name); 
           Kgeo_680nm=double(Kgeo_680nm')*BRDF_factor;
           Dataset_name=['/',tm{hi},'/Kiso_680nm'];
           Kiso_680nm=h5read([filepath,filestr(fi).name],Dataset_name); 
           Kiso_680nm=double(Kiso_680nm')*BRDF_factor;
           Dataset_name=['/',tm{hi},'/Kvol_680nm'];
           Kvol_680nm=h5read([filepath,filestr(fi).name],Dataset_name); 
           Kvol_680nm=double(Kvol_680nm')*BRDF_factor;
           
              Dataset_name=['/',tm{hi},'/Kgeo_780nm'];
           Kgeo_780nm=h5read([filepath,filestr(fi).name],Dataset_name); 
           Kgeo_780nm=double(Kgeo_780nm')*BRDF_factor;
           Dataset_name=['/',tm{hi},'/Kiso_780nm'];
           Kiso_780nm=h5read([filepath,filestr(fi).name],Dataset_name); 
           Kiso_780nm=double(Kiso_780nm')*BRDF_factor;
           Dataset_name=['/',tm{hi},'/Kvol_780nm'];
           Kvol_780nm=h5read([filepath,filestr(fi).name],Dataset_name); 
           Kvol_780nm=double(Kvol_780nm')*BRDF_factor;

           Dataset_name=['/',tm{hi},'/LW_mask'];
           LW_mask=h5read([filepath,filestr(fi).name],Dataset_name); 
           LW_mask  = LW_mask';
          
              Dataset_name=['/',tm{hi},'/Status_QA'];
           Status_QA=h5read([filepath,filestr(fi).name],Dataset_name); 
           Status_QA  = Status_QA';


          % id=find(Red<=0 | NIR<=0 | VZA1<0 | SZA1<0 | VAA1<0 | SAA1<0);
          % Red(id)=NaN;
          % NIR(id)=NaN;
          % VZA1(id)=NaN;
          % SZA1(id)=NaN;
          % VAA1(id)=NaN;
          % SAA1(id)=NaN;

              % Output  
   
    fileoutname=[fileoutpath,filename_tmp(1:end-6),'_Refl.mat'];
   save(fileoutname,'Red','NIR','VZA1','SZA1','VAA1','SAA1','Kgeo_443nm','Kiso_443nm','Kvol_443nm',...
     'Kgeo_551nm','Kiso_551nm','Kvol_551nm','Kgeo_680nm','Kiso_680nm','Kvol_680nm','Kgeo_780nm','Kiso_780nm','Kvol_780nm',...
        'AOD','Cloud_mask','LW_mask','Status_QA');

   disp(['Month:',num2str(mi),'   Day:',num2str(dayi)]);
        end        
    end
   
    end
%end
end
