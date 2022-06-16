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

%%

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

R_Cutoff_Value = 11000;
G_Cutoff_Value = 6000;
B_Cutoff_Value = 15000;

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
%%
figure(5)
subplot(1,2,1)
imshow(R_Tiff_Xnm_Crop) 
subplot(1,2,2)
imshow(R_Tiff_0nm_Crop)

figure(6)
subplot(1,2,1)
imshow(G_Tiff_Xnm_Crop) 
subplot(1,2,2)
imshow(G_Tiff_0nm_Crop)

figure(7)
subplot(1,2,1)
imshow(B_Tiff_Xnm_Crop) 
subplot(1,2,2)
imshow(B_Tiff_0nm_Crop)
%%

% mean_R_0nm = mean(mean(nonzeros(R_Tiff_0nm_Crop)));
% mean_G_0nm = mean(mean(nonzeros(G_Tiff_0nm_Crop)));
% mean_B_0nm = mean(mean(nonzeros(B_Tiff_0nm_Crop)));
% mean_R_Xnm = mean(mean(nonzeros(R_Tiff_Xnm_Crop)));
% mean_G_Xnm = mean(mean(nonzeros(G_Tiff_Xnm_Crop)));
% mean_B_Xnm = mean(mean(nonzeros(B_Tiff_Xnm_Crop)));

% I_inR = (double(mean_R_0nm) ./ Ref_at_0nm_R);
% I_inG = (double(mean_G_0nm) ./ Ref_at_0nm_G);
% I_inB = (double(mean_B_0nm) ./ Ref_at_0nm_B);
% 
% RefRed_at_XnmO2 = (double(mean_R_Xnm)./I_inR);
% RefGre_at_XnmO2 = (double(mean_G_Xnm)./I_inG);
% RefBlu_at_XnmO2 = (double(mean_B_Xnm)./I_inB);

numel_R_Tiff_Xnm_Crop = sum(sum(R_Tiff_Xnm_Crop>0));
numel_G_Tiff_Xnm_Crop = sum(sum(G_Tiff_Xnm_Crop>0));
numel_B_Tiff_Xnm_Crop = sum(sum(B_Tiff_Xnm_Crop>0));

numel_R_Tiff_0nm_Crop = sum(sum(R_Tiff_0nm_Crop>0));
numel_G_Tiff_0nm_Crop = sum(sum(G_Tiff_0nm_Crop>0));
numel_B_Tiff_0nm_Crop = sum(sum(B_Tiff_0nm_Crop>0));
%
[rowsR, colsR, ~] = size(R_Tiff_Bothnm_Crop);
[rowsG, colsG, ~] = size(G_Tiff_Bothnm_Crop);
[rowsB, colsB, ~] = size(B_Tiff_Bothnm_Crop);
%
R_Center_x = round((1+rowsR)/2);
R_Center_y = round((1+colsR)/2);
G_Center_x = round((1+rowsG)/2);
G_Center_y = round((1+colsG)/2);
B_Center_x = round((1+rowsB)/2);
B_Center_y = round((1+colsB)/2);

numel_crop_R = numel_R_Tiff_0nm_Crop;
numel_crop_G = numel_G_Tiff_0nm_Crop;
numel_crop_B = numel_B_Tiff_0nm_Crop;

rect3_R_w = round(sqrt(numel_crop_R));
rect4_R_h = round(sqrt(numel_crop_R));
rect3_G_w = round(sqrt(numel_crop_G));
rect4_G_h = round(sqrt(numel_crop_G));
rect3_B_w = round(sqrt(numel_crop_B));
rect4_B_h = round(sqrt(numel_crop_B));

rect_R = [R_Center_x-50 R_Center_y-50  rect3_R_w rect4_R_h];
rect_G = [G_Center_x-50 G_Center_y-50  rect3_G_w rect4_G_h]; %DO NOT HARD CODE CHANGE -50
rect_B = [B_Center_x-50 B_Center_y-50  rect3_B_w rect4_B_h];

cropped_im_R_Xnm = imcrop(R_Tiff_Xnm_Crop, rect_R);
cropped_im_G_Xnm = imcrop(G_Tiff_Xnm_Crop, rect_G);
cropped_im_B_Xnm = imcrop(B_Tiff_Xnm_Crop, rect_B);

figure(8)
imshow(cropped_im_R_Xnm)
title('Cropped R Xnm ROI')
figure(9)
imshow(cropped_im_G_Xnm)
title('Cropped G Xnm ROI')
figure(10)
imshow(cropped_im_B_Xnm)
title('Cropped B Xnm ROI')

numel_crop_R_real = numel(cropped_im_R_Xnm)
numel_crop_G_real = numel(cropped_im_G_Xnm)
numel_crop_B_real = numel(cropped_im_B_Xnm)

nonzeros_R_0nm = nonzeros(R_Tiff_0nm_Crop); %numel = 15129
nonzeros_R_Xnm = nonzeros(cropped_im_R_Xnm); %numel = 14851
nonzeros_G_0nm = nonzeros(G_Tiff_0nm_Crop); %numel = 15129
nonzeros_G_Xnm = nonzeros(cropped_im_G_Xnm); %numel = 14851
nonzeros_B_0nm = nonzeros(B_Tiff_0nm_Crop); %numel = 15129
nonzeros_B_Xnm = nonzeros(cropped_im_B_Xnm); %numel = 14851

if numel(nonzeros_R_0nm) > numel(nonzeros_R_Xnm)
    nonzeros_R_0nm = nonzeros_R_0nm(1:numel(nonzeros_R_Xnm))
elseif numel(nonzeros_R_0nm) < numel(nonzeros_R_Xnm)
    nonzeros_R_Xnm = nonzeros_R_Xnm(1:numel(nonzeros_R_0nm)) 
end

if numel(nonzeros_G_0nm) > numel(nonzeros_G_Xnm)
    nonzeros_G_0nm = nonzeros_G_0nm(1:numel(nonzeros_G_Xnm))
elseif numel(nonzeros_G_0nm) < numel(nonzeros_G_Xnm)
    nonzeros_G_Xnm = nonzeros_G_Xnm(1:numel(nonzeros_G_0nm)) 
end

if numel(nonzeros_B_0nm) > numel(nonzeros_B_Xnm)
    nonzeros_B_0nm = nonzeros_B_0nm(1:numel(nonzeros_B_Xnm))
elseif numel(nonzeros_B_0nm) < numel(nonzeros_B_Xnm)
    nonzeros_B_Xnm = nonzeros_B_Xnm(1:numel(nonzeros_B_0nm)) 
end

%%


% [rows, columns, ~] = size(R_Tiff_Bothnm_Crop)
% R_Tiff_Xnm_Crop = zeros(rows,columns) 

% [R_Tiff_Xnm_Crop,x_crop_Xnm,y_crop_Xnm] = AreaSelection_Circle(R_Tiff);
% G_Tiff_Xnm_Crop = AreaSelection_Circle_Mod(G_Tiff,x_crop_Xnm,y_crop_Xnm);
% B_Tiff_Xnm_Crop = AreaSelection_Circle_Mod(B_Tiff,x_crop_Xnm,y_crop_Xnm);

% R_Tiff_0nm_Crop = R_Tiff_Bothnm_Crop - R_Tiff_Xnm_Crop



%Intensity In calculation (calibration) I_in = I_out_Si / Reflectance
%
I_inR = round(double(nonzeros_R_0nm) ./ Ref_at_0nm_R);
I_inG = round(double(nonzeros_G_0nm) ./ Ref_at_0nm_G);
I_inB = round(double(nonzeros_B_0nm) ./ Ref_at_0nm_B);
%%
%Reflectance from oxide calculation: Reflectance_Oxide = I_Out_SiO2 / I_in
%Multiplied with masks to eliminate other pixels and thus divide by zeros.

RefRed_at_XnmO2 = (double(nonzeros_R_Xnm)./I_inR);
RefGre_at_XnmO2 = (double(nonzeros_G_Xnm)./I_inG);
RefBlu_at_XnmO2 = (double(nonzeros_B_Xnm)./I_inB);
%%
%Take the average of Reflectance values in R,G,B only pixels

MeanRefRed_at_XnmO2 = mean(mean(RefRed_at_XnmO2));
MeanRefGre_at_XnmO2 = mean(mean(RefGre_at_XnmO2));
MeanRefBlu_at_XnmO2 = mean(mean(RefBlu_at_XnmO2));

%use reftocurve function to estimate thickness from these reflectance values
esti_L = reftocurve(MeanRefRed_at_XnmO2,MeanRefGre_at_XnmO2,MeanRefBlu_at_XnmO2);
esti_L_lsqr = reftocurve_lsqr(MeanRefRed_at_XnmO2,MeanRefGre_at_XnmO2,MeanRefBlu_at_XnmO2);

% esti_L = reftocurve(RefRed_at_XnmO2,RefGre_at_XnmO2,RefBlu_at_XnmO2);
% esti_L_lsqr = reftocurve_lsqr(RefRed_at_XnmO2,RefGre_at_XnmO2,RefBlu_at_XnmO2);

figure(15)

hold on
plot(lambda,Gamma(:,valtoindex_L(abs(esti_L_lsqr))),'m--','LineWidth',2)           %Estimated nm thickness
plot(lambda,Gamma(:,valtoindex_L(0.12)),'c--','LineWidth',2)                       %120 nm thickness

plot(cw_r,MeanRefRed_at_XnmO2,'r.','MarkerSize',30)
plot(cw_g,MeanRefGre_at_XnmO2,'g.','MarkerSize',30)
plot(cw_b,MeanRefBlu_at_XnmO2,'b.','MarkerSize',30)

% plot(cw_r,RefRed_at_XnmO2,'r.','MarkerSize',30)
% plot(cw_g,RefGre_at_XnmO2,'g.','MarkerSize',30)
% plot(cw_b,RefBlu_at_XnmO2,'b.','MarkerSize',30)

xlabel('lambda (\mum)')
ylabel('Reflectance')
xlim([0.4 0.68])
ylim([0 1])
legend1 = sprintf('L = %0.2f nm', 1000*esti_L_lsqr);
legend2 = sprintf('L = 120 nm');
legend(legend1, legend2,'Ref at R','Ref at G','Ref at B','location','bestoutside')

if isnumeric(theoric_L)
percent_error = 100*((abs(esti_L_lsqr-theoric_L))/theoric_L);
fprintf('Actual thickness is %f nm \n',theoric_L*1000)
fprintf('Estimated thickness is %f nm \n',esti_L_lsqr*1000)
fprintf("Percent error is %f", percent_error)
elseif theoric_L == 'UNKNOWN'
fprintf('Estimated thickness is %f nm \n',esti_L_lsqr*1000)
disp('Percent Error is not known')
end
% 
% %branch test








