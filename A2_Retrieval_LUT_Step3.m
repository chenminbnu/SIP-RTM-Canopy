function [retrieval_value] = A2_Retrieval_LUT_Step3(Refs_red, Refs_nir,SZAs, VZAs, RAAs,RSO_LUT)

u = size(Refs_red,1);
INPUT=zeros(u,8);
INPUT(:,1)=SZAs;  %SZA
INPUT(:,2)=VZAs;  %VZA
INPUT(:,3)=RAAs;  %VZA

INPUT(:,4)=0.5*ones(u,1);
INPUT(:,5)=0.5*ones(u,1);
INPUT(:,6)=0.1*ones(u,1);
INPUT(:,7)=1*ones(u,1);
INPUT(:,8)=1*ones(u,1);
[location1]=findlocation_RSO(INPUT);
nlength=6*6*70*8*10;
errors=zeros(u,nlength);
delta2_red = 0.005;
delta2_NIR = 0.014;
for i2=1:nlength
    for i3=1:u
    %%   NIRv_epic = (Refs_nir(i3,1)-Refs_red(i3,1))./(Refs_nir(i3,1)+Refs_red(i3,1)).*Refs_nir(i3,1);
      %    NIRv_lut = (RSO_LUT(location1(i3,1)+i2-1,2)-RSO_LUT(location1(i3,1)+i2-1,1))./(RSO_LUT(location1(i3,1)+i2-1,2)+RSO_LUT(location1(i3,1)+i2-1,1)).*RSO_LUT(location1(i3,1)+i2-1,2);

          errors(i3,i2)=(RSO_LUT(location1(i3,1)+i2-1,1)-Refs_red(i3,1))^2/delta2_red+(RSO_LUT(location1(i3,1)+i2-1,2)-Refs_nir(i3,1))^2/delta2_NIR;
         %% errors(i3,i2)=( NIRv_lut-NIRv_epic)^2/delta2_NIR;

    end
end
errors = sum(errors,1);
errors = errors';
[error_rows]=find(errors<2*u);
%[error_rows]=find(errors<1*u);

retrieval_value= zeros(5,1);
if(isempty(error_rows))
;
else
error_values = 1./(errors(error_rows)+0.000001);
%CIy1 CIy2 LAI LIDFa FPAR i0 NIRv

[OUTPUT]=findvalue_RSO(error_rows);



for i = 1:5
    retrieval_value(i,1) = sum(OUTPUT(:,i).*(error_values))/sum(error_values);
end
end
end






