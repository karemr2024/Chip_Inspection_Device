clc; close all; clearvars;

%Below are prompts for user to enter unknown SiO2 nm (120 nm for now)
%and 0 nm (Si) chip images.

[R_Tiff_Name,R_Tiff_Path] = uigetfile('*.tif','Red Image'); %Import Unknown Thickness: Red Image
if isequal(R_Tiff_Name,0)
   disp('User selected Cancel');
else
   disp(['User selected ', fullfile(R_Tiff_Name,R_Tiff_Path)]);
   R_Tiff = imread(strcat(R_Tiff_Path, R_Tiff_Name));
end

[G_Tiff_Name,G_Tiff_Path] = uigetfile('*.tif','Green Image'); %Import Unknown Thickness: Green Image
if isequal(G_Tiff_Name,0)
   disp('User selected Cancel');
else
   disp(['User selected ', fullfile(G_Tiff_Name,G_Tiff_Path)]);
   G_Tiff = imread(strcat(G_Tiff_Path, G_Tiff_Name));
end

[B_Tiff_Name,B_Tiff_Path] = uigetfile('*.tif','Blue Image'); %Import Unknown Thickness: Blue Image 
if isequal(B_Tiff_Name,0)
   disp('User selected Cancel');
else
   disp(['User selected ', fullfile(B_Tiff_Name,B_Tiff_Path)]);
   B_Tiff = imread(strcat(B_Tiff_Path, B_Tiff_Name));
end

%%
clearvars -except B_Tiff G_Tiff R_Tiff

prompt_NorS = input('Do you know the thickness of your chip? (Y/N)\n', "s");
if prompt_NorS == 'Y'
theoric_L = (input('Enter the thickness of your chip in nm: \n'))/1000;
else
theoric_L = 'UNKNOWN';
end

load("Imaging_Data.mat")

RGB_Combo = R_Tiff + G_Tiff + B_Tiff; %Add R,G,B Channels to create a RGB combo image

%Reflectances at 0 and 120 nm

%ADD THESE TO IMAGING_DATA.MAT
Ref_at_0nm_R = Gamma(valtoindex_lambda(cw_r),valtoindex_L(0));
Ref_at_0nm_G = Gamma(valtoindex_lambda(cw_g),valtoindex_L(0));
Ref_at_0nm_B = Gamma(valtoindex_lambda(cw_b),valtoindex_L(0));
Ref_at_120nm_R = Gamma(valtoindex_lambda(cw_r),valtoindex_L(0.12));
Ref_at_120nm_G = Gamma(valtoindex_lambda(cw_g),valtoindex_L(0.12));
Ref_at_120nm_B = Gamma(valtoindex_lambda(cw_b),valtoindex_L(0.12));

%Jun 14 Note: Can't implement AreaSelection_Circle function!! 

[R_Tiff_Bothnm_Crop,x_crop_Bothnm,y_crop_Bothnm] = AreaSelection_Circle(R_Tiff);
G_Tiff_Bothnm_Crop = AreaSelection_Circle_Mod(G_Tiff,x_crop_Bothnm,y_crop_Bothnm);
B_Tiff_Bothnm_Crop = AreaSelection_Circle_Mod(B_Tiff,x_crop_Bothnm,y_crop_Bothnm);

R_Cutoff_Value = 16000;
G_Cutoff_Value = 6000;
B_Cutoff_Value = 8000;

R_Cutoff_Mask_0nm = uint16(R_Tiff_Bothnm_Crop > R_Cutoff_Value);
G_Cutoff_Mask_0nm = uint16(G_Tiff_Bothnm_Crop > G_Cutoff_Value);
B_Cutoff_Mask_0nm = uint16(B_Tiff_Bothnm_Crop > B_Cutoff_Value);

R_Cutoff_Mask_Xnm = uint16(R_Tiff_Bothnm_Crop < R_Cutoff_Value);
G_Cutoff_Mask_Xnm = uint16(G_Tiff_Bothnm_Crop < G_Cutoff_Value);
B_Cutoff_Mask_Xnm = uint16(B_Tiff_Bothnm_Crop < B_Cutoff_Value);

R_Tiff_0nm_Crop = R_Tiff_Bothnm_Crop.*R_Cutoff_Mask_0nm;
G_Tiff_0nm_Crop = G_Tiff_Bothnm_Crop.*G_Cutoff_Mask_0nm;
B_Tiff_0nm_Crop = B_Tiff_Bothnm_Crop.*B_Cutoff_Mask_0nm;

R_Tiff_Xnm_Crop = R_Tiff_Bothnm_Crop.*R_Cutoff_Mask_Xnm;
G_Tiff_Xnm_Crop = G_Tiff_Bothnm_Crop.*G_Cutoff_Mask_Xnm;
B_Tiff_Xnm_Crop = B_Tiff_Bothnm_Crop.*B_Cutoff_Mask_Xnm;
%
R_Tiff_0nm_Crop_nonzeros = R_Tiff_0nm_Crop > 0;
numel_R_Tiff_Xnm_Crop = sum(sum(R_Cutoff_Mask_Xnm))
numel_G_Tiff_Xnm_Crop = sum(sum(G_Cutoff_Mask_Xnm))
numel_B_Tiff_Xnm_Crop = sum(sum(B_Cutoff_Mask_Xnm))
%
[rows, cols, ~] = size(R_Tiff_Bothnm_Crop);
%
R_Center_x = round((1+rows)/2)
R_Center_y = round((1+cols)/2)

rect = [R_Center_x R_Center_x  64 64]
%

subplot(1,2,1)
imshow(R_Tiff_Xnm_Crop)
subplot(1,2,2)
imshow(R_Tiff_0nm_Crop)
% [rows, columns, ~] = size(R_Tiff_Bothnm_Crop)
% R_Tiff_Xnm_Crop = zeros(rows,columns) 

% [R_Tiff_Xnm_Crop,x_crop_Xnm,y_crop_Xnm] = AreaSelection_Circle(R_Tiff);
% G_Tiff_Xnm_Crop = AreaSelection_Circle_Mod(G_Tiff,x_crop_Xnm,y_crop_Xnm);
% B_Tiff_Xnm_Crop = AreaSelection_Circle_Mod(B_Tiff,x_crop_Xnm,y_crop_Xnm);

% R_Tiff_0nm_Crop = R_Tiff_Bothnm_Crop - R_Tiff_Xnm_Crop



%Intensity In calculation (calibration) I_in = I_out_Si / Reflectance
%%
I_inR = round(double(R_Tiff_0nm_Crop) ./ Ref_at_0nm_R);
I_inG = round(double(G_Tiff_0nm_Crop) ./ Ref_at_0nm_G);
I_inB = round(double(B_Tiff_0nm_Crop) ./ Ref_at_0nm_B);
%%
%Reflectance from oxide calculation: Reflectance_Oxide = I_Out_SiO2 / I_in
%Multiplied with masks to eliminate other pixels and thus divide by zeros.

RefRed_at_XnmO2 = (double(R_Tiff_Xnm_Crop)./I_inR);
RefGre_at_XnmO2 = (double(G_Tiff_Xnm_Crop)./I_inG);
RefBlu_at_XnmO2 = (double(B_Tiff_Xnm_Crop)./I_inB);
%%
% %Take the average of Reflectance values in R,G,B only pixels

MeanRefRed_at_XnmO2 = mean(mean(RefRed_at_XnmO2));
MeanRefGre_at_XnmO2 = mean(mean(RefGre_at_XnmO2));
MeanRefBlu_at_XnmO2 = mean(mean(RefBlu_at_XnmO2));

%use reftocurve function to estimate thickness from these reflectance values
esti_L = reftocurve(MeanRefRed_at_XnmO2,MeanRefGre_at_XnmO2,MeanRefBlu_at_XnmO2);
esti_L_lsqr = reftocurve_lsqr(MeanRefRed_at_XnmO2,MeanRefGre_at_XnmO2,MeanRefBlu_at_XnmO2);

% figure(15)
% 
% hold on
% plot(lambda,Gamma(:,valtoindex_L(abs(esti_L_lsqr))),'m--','LineWidth',2)           %Estimated nm thickness
% plot(lambda,Gamma(:,valtoindex_L(0.12)),'c--','LineWidth',2)                       %120 nm thickness
% 
% plot(cw_r,MeanRefRed_at_XnmO2,'r.','MarkerSize',30)
% plot(cw_g,MeanRefGre_at_XnmO2,'g.','MarkerSize',30)
% plot(cw_b,MeanRefBlu_at_XnmO2,'b.','MarkerSize',30)
% 
% xlabel('lambda (\mum)')
% ylabel('Reflectance')
% xlim([0.4 0.68])
% ylim([0 1])
% legend1 = sprintf('L = %0.2f nm', 1000*esti_L_lsqr);
% legend2 = sprintf('L = 120 nm');
% legend(legend1, legend2,'Ref at R','Ref at G','Ref at B','location','bestoutside')
% 
% if isnumeric(theoric_L)
% percent_error = 100*((abs(esti_L_lsqr-theoric_L))/theoric_L);
% fprintf('Actual thickness is %f nm \n',theoric_L*1000)
% fprintf('Estimated thickness is %f nm \n',esti_L_lsqr*1000)
% fprintf("Percent error is %f", percent_error)
% elseif theoric_L == 'UNKNOWN'
% fprintf('Estimated thickness is %f nm \n',esti_L_lsqr*1000)
% disp('Percent Error is not known')
% end
% 
% %branch test







