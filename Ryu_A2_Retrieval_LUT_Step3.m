% %Hemiview
% clear
% % RD_LUT=zeros(11*7*6*6*70*8,8);
% load RD_LUT.mat
load INPUT_DATA.mat % time diffuseratio SZA NIRv LAI FPAR

[u0,v0]=size(INPUT_DATA);
OUTPUT_DATA=zeros(u0,14);

for i1=132:244
    
var_temp={find(INPUT_DATA(:,1)>max(i1-2,132)),find(INPUT_DATA(:,1)<min(i1+3,245))};   %5 days
ind_filter=intersectvecs(var_temp);
[u,v]=size(ind_filter);
INPUT_temp=INPUT_DATA(ind_filter,:);

var_temp2={find(INPUT_DATA(:,1)>i1),find(INPUT_DATA(:,1)<i1+1)}; 
ind_filter2=intersectvecs(var_temp2);
[u2,v2]=size(ind_filter2);

[C,ia,ib] = intersect(ind_filter,ind_filter2);

    INPUT=zeros(u,1);
    INPUT(:,1)=INPUT_temp(:,2);  %diffuse ratio
    INPUT(:,2)=INPUT_temp(:,3);  %SZA

    INPUT(:,3)=0.5*ones(u,1);
    INPUT(:,4)=0.5*ones(u,1);
    INPUT(:,5)=0.1*ones(u,1);
    INPUT(:,6)=10*ones(u,1);
    [location1]=findlocation_RD(INPUT);
    
    nlength=6*6*70*8;
    error=zeros(u,nlength);
    num=10;
    
    for i2=1:nlength
        for i3=1:u
        error(i3,i2)=abs(RD_LUT(location1(i3,1)+i2-1,8)-INPUT_temp(i3,4));
        end
    end
    
    error_med=median(error);
    error_med=error_med';
    [error_value, error_row]=sort(error_med);
        
    %CIy1 CIy2 LAI LIDFa FPAR i0 NIRv    
    
    [OUTPUT]=findvalue(error_row(1:num,1));
    
    for j1=1:4
        OUTPUT_DATA(ind_filter2,j1)=mean(OUTPUT(:,j1));
        OUTPUT_DATA(ind_filter2,j1+7)=std(OUTPUT(:,j1));
    end
    
    for j2=1:3
        for j3=1:u2
        OUTPUT_DATA(ind_filter(ia(j3)),j2+4)=mean(RD_LUT(location1(ia(j3))+error_row(1:num,1)-1,j2+5));
        OUTPUT_DATA(ind_filter(ia(j3)),j2+4+7)=std(RD_LUT(location1(ia(j3))+error_row(1:num,1)-1,j2+5));
        end
    end
end

    % load INPUT_DATA.mat % time diffuseratio SZA NIRv LAI FPAR
%CIy1 CIy2 LAI LIDFa FPAR i0 NIRv


subplot (3,2,1)
plot(INPUT_DATA(:,1),INPUT_DATA(:,4),'r.',INPUT_DATA(:,1),INPUT_DATA(:,6),'b.',INPUT_DATA(:,1),OUTPUT_DATA(:,5),'g.');
xlabel('DOY','fontsize',13);
ylabel('Value','fontsize',13);
legend('Field NIRv','Field FPAR','Retrieved FPAR');

subplot (3,2,2)
plot(INPUT_DATA(:,1),OUTPUT_DATA(:,3),'r*',INPUT_DATA(:,1),INPUT_DATA(:,5),'b*');
xlabel('DOY','fontsize',13);
ylabel('Value','fontsize',13);
legend('Retrieved LAI','Field LAI');

subplot (3,2,3)
plot(INPUT_DATA(:,1),OUTPUT_DATA(:,1),'b*');  %INPUT_DATA(:,1),OUTPUT_DATA(:,2),'b.'
xlabel('DOY','fontsize',13);
ylabel('Clumping index (CI)','fontsize',13);
legend('Retrieved CI at 20 degree');

subplot (3,2,4)
plot(INPUT_DATA(:,1),OUTPUT_DATA(:,4),'b*');
xlabel('DOY','fontsize',13);
ylabel('Mean leaf angle','fontsize',13);
legend('Retrieved mean leaf angle');

% subplot (3,2,5)
% plot(INPUT_DATA(:,4),INPUT_DATA(:,6),'b.',OUTPUT_DATA(:,7),OUTPUT_DATA(:,5),'r.');
% xlabel('NIRv','fontsize',13);
% ylabel('FPAR','fontsize',13);
% legend('Field','Retrieved');

subplot (3,2,5)
plot(OUTPUT_DATA(:,6),OUTPUT_DATA(:,5),'b*');
xlabel('Retrieved i0','fontsize',13);
ylabel('Retrieved FPAR','fontsize',13);
legend('FPAR vs. i0');


subplot (3,2,6)
plot(INPUT_DATA(:,6),OUTPUT_DATA(:,5),'b.');
xlabel('Field FPAR','fontsize',13);
ylabel('Retrieved FPAR','fontsize',13);
% legend('NIRv','FPAR');


title('Preliminary performance at the rice paddy site');








% Directional view
% % RO_LUT=zeros(11*7*7*7*6*6*70*8,8);
% 
% data_out=zeros(npixel,12);
% for i1=1:npixel
%     INPUT(1)=Ratio;
%     INPUT(2)=tts;
%     INPUT(3)=tto;
%     INPUT(4)=psi;
%     INPUT(4)=abs(psi-360*round(psi/360));  %           (to ensure that volscatt is symmetric for psi=90 and psi=270)
%     
%     INPUT(5)=0.5;
%     INPUT(6)=0.5;
%     INPUT(7)=0.1;
%     INPUT(8)=10;
%     [location1]=findlocation_RO(INPUT);
%     
%     nlength=6*6*70*8;
%     error=zeros(nlength);
%     num=10;
%     
%     for i2=1:nlength
%         error(i2)=abs(RO_LUT(location1+i2-1,8)-NIRv);
%     end
%     [error_value, error_row]=sort(error(:));
%     
%     
%     %CIy1 CIy2 LAI LIDFa FPAR i0 NIRv
%     
%     [OUTPUT]=findvalue_RO(error_row);
%     
%     data_out(i1,1)=mean(OUTPUT(:,1));
%     data_out(i1,2)=std(OUTPUT(:,1));
%     data_out(i1,3)=mean(OUTPUT(:,2));
%     data_out(i1,4)=std(OUTPUT(:,2));
%     data_out(i1,5)=mean(OUTPUT(:,3));
%     data_out(i1,6)=std(OUTPUT(:,3));
%     data_out(i1,7)=mean(OUTPUT(:,4));
%     data_out(i1,8)=std(OUTPUT(:,4));
%     data_out(i1,9)=mean(RO_LUT(location1+error_row-1,6));
%     data_out(i1,10)=std(RO_LUT(location1+error_row-1,6));
%     data_out(i1,11)=mean(RO_LUT(location1+error_row-1,7));
%     data_out(i1,12)=std(RO_LUT(location1+error_row-1,7));
% end




