function [LAI] = RETC_Retrieval(Refs_red, Refs_nir,RSO_LUT)

nlength=70;
errors=zeros(nlength,1);
errors_red=zeros(nlength,1);
errors_NIR=zeros(nlength,1);

delta2_red = 0.005;
delta2_NIR = 0.014;
for i2=1:nlength
    %%   NIRv_epic = (Refs_nir(i3,1)-Refs_red(i3,1))./(Refs_nir(i3,1)+Refs_red(i3,1)).*Refs_nir(i3,1);
      %    NIRv_lut = (RSO_LUT(location1(i3,1)+i2-1,2)-RSO_LUT(location1(i3,1)+i2-1,1))./(RSO_LUT(location1(i3,1)+i2-1,2)+RSO_LUT(location1(i3,1)+i2-1,1)).*RSO_LUT(location1(i3,1)+i2-1,2);
          errors(i2,1)=(RSO_LUT(i2,1)-Refs_red)^2/delta2_red+(RSO_LUT(i2,2)-Refs_nir)^2/delta2_NIR;
              errors_red(i2,1)=RSO_LUT(i2,1)-Refs_red;
              errors_NIR(i2,1)=RSO_LUT(i2,2)-Refs_nir;

          %% errors(i3,i2)=( NIRv_lut-NIRv_epic)^2/delta2_NIR;
end
[error_rows]=find(errors<2);
%[error_rows]=find(errors<1*u);
%   scatter(errors_red,errors_NIR)
%   title('forest')
if(isempty(error_rows))
LAI=0;
else
error_values = 1./(errors(error_rows)+0.000001);
%CIy1 CIy2 LAI LIDFa FPAR i0 NIRv

LAIs=error_rows*0.1;
  LAI = sum(LAIs.*error_values)/sum(error_values);
end
end






