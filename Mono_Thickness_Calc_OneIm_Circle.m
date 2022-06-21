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

[R_Tiff_Bothnm_Crop,x_crop_Bothnm,y_crop_Bothnm] = AreaSelection_Circle(R_Tiff);
G_Tiff_Bothnm_Crop = AreaSelection_Circle_Mod(G_Tiff,x_crop_Bothnm,y_crop_Bothnm);
B_Tiff_Bothnm_Crop = AreaSelection_Circle_Mod(B_Tiff,x_crop_Bothnm,y_crop_Bothnm);

[cropped_im_R_0nm_out,x_crop_0nm_out,y_crop_0nm_out] = AreaSelection_Circle_SiOut(R_Tiff_Bothnm_Crop);
cropped_im_G_0nm_out = AreaSelection_Circle_Mod(G_Tiff_Bothnm_Crop,x_crop_0nm_out,y_crop_0nm_out);
cropped_im_B_0nm_out = AreaSelection_Circle_Mod(B_Tiff_Bothnm_Crop,x_crop_0nm_out,y_crop_0nm_out);

[cropped_im_R_0nm_in,x_crop_0nm_in,y_crop_0nm_in] = AreaSelection_Circle_SiIn(R_Tiff_Bothnm_Crop,x_crop_0nm_out,y_crop_0nm_out);
cropped_im_G_0nm_in = AreaSelection_Circle_Mod(G_Tiff_Bothnm_Crop,x_crop_0nm_in,y_crop_0nm_in);
cropped_im_B_0nm_in = AreaSelection_Circle_Mod(B_Tiff_Bothnm_Crop,x_crop_0nm_in,y_crop_0nm_in);

[rows_0nm_out,cols_0nm_out] = size(cropped_im_R_0nm_out)
[rows_0nm_in,cols_0nm_in] = size(cropped_im_R_0nm_in)
rowdiff = rows_0nm_out(1) - rows_0nm_in(1);
coldiff = cols_0nm_out(1) - cols_0nm_in(1);
threshold_rows = round((rowdiff/2)+1)
threshold_cols = round((coldiff/2)+1)
%%
bigger_120nm_R(threshold:rows_0nm_in+threshold,threshold:rows_0nm_in+threshold) = cropped_im_R_0nm_in(1:rows_0nm_in+1,1:rows_0nm_in+1)   
%FIX HERE
%%

%%Different sizes! try automatically drawing the inner circle after drawing
%%outer circle.



% cropped_Ring_R = 


%%

% R_Cutoff_Value = (max(max(cropped_im_R_Xnm))+max(max(R_Tiff_Bothnm_Crop)))/2;
% G_Cutoff_Value = (max(max(cropped_im_G_Xnm))+max(max(G_Tiff_Bothnm_Crop)))/2;
% B_Cutoff_Value = (max(max(cropped_im_B_Xnm))+max(max(B_Tiff_Bothnm_Crop)))/2;

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



numel_R_Tiff_Xnm_Crop = sum(sum(R_Tiff_Xnm_Crop>0));
numel_G_Tiff_Xnm_Crop = sum(sum(G_Tiff_Xnm_Crop>0));
numel_B_Tiff_Xnm_Crop = sum(sum(B_Tiff_Xnm_Crop>0));

numel_R_Tiff_0nm_Crop = sum(sum(R_Tiff_0nm_Crop>0));
numel_G_Tiff_0nm_Crop = sum(sum(G_Tiff_0nm_Crop>0));
numel_B_Tiff_0nm_Crop = sum(sum(B_Tiff_0nm_Crop>0));

[rowsR, colsR, ~] = size(R_Tiff_Bothnm_Crop);
[rowsG, colsG, ~] = size(G_Tiff_Bothnm_Crop);
[rowsB, colsB, ~] = size(B_Tiff_Bothnm_Crop);

nonzeros_R_0nm = nonzeros(R_Tiff_0nm_Crop); 
nonzeros_R_Xnm = nonzeros(cropped_im_R_Xnm); 
nonzeros_G_0nm = nonzeros(G_Tiff_0nm_Crop); 
nonzeros_G_Xnm = nonzeros(cropped_im_G_Xnm); 
nonzeros_B_0nm = nonzeros(B_Tiff_0nm_Crop);
nonzeros_B_Xnm = nonzeros(cropped_im_B_Xnm); 

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

%Intensity In calculation (calibration) I_in = I_out_Si / Reflectance
%
I_inR = round(double(nonzeros_R_0nm) ./ Ref_at_0nm_R);
I_inG = round(double(nonzeros_G_0nm) ./ Ref_at_0nm_G);
I_inB = round(double(nonzeros_B_0nm) ./ Ref_at_0nm_B);


%Reflectance from oxide calculation: Reflectance_Oxide = I_Out_SiO2 / I_in
%Multiplied with masks to eliminate other pixels and thus divide by zeros.

RefRed_at_XnmO2 = (double(nonzeros_R_Xnm)./I_inR);
RefGre_at_XnmO2 = (double(nonzeros_G_Xnm)./I_inG);
RefBlu_at_XnmO2 = (double(nonzeros_B_Xnm)./I_inB);
%
%Take the average of Reflectance values in R,G,B only pixels

MeanRefRed_at_XnmO2 = mean(mean(RefRed_at_XnmO2));
MeanRefGre_at_XnmO2 = mean(mean(RefGre_at_XnmO2));
MeanRefBlu_at_XnmO2 = mean(mean(RefBlu_at_XnmO2));

%use reftocurve function to estimate thickness from these reflectance values
esti_L = reftocurve(MeanRefRed_at_XnmO2,MeanRefGre_at_XnmO2,MeanRefBlu_at_XnmO2);
esti_L_lsqr = reftocurve_lsqr(MeanRefRed_at_XnmO2,MeanRefGre_at_XnmO2,MeanRefBlu_at_XnmO2);
%%
figure(15)

%ADD TITLE AS INPUT FROM USER

hold on
plot(lambda,Gamma(:,valtoindex_L(abs(esti_L_lsqr))),'m--','LineWidth',2)           %Estimated nm thickness
plot(lambda,Gamma(:,valtoindex_L(theoric_L)),'c--','LineWidth',2)                       %120 nm thickness

plot(cw_r,MeanRefRed_at_XnmO2,'r.','MarkerSize',30)
plot(cw_g,MeanRefGre_at_XnmO2,'g.','MarkerSize',30)
plot(cw_b,MeanRefBlu_at_XnmO2,'b.','MarkerSize',30)

xlabel('lambda (\mum)')
ylabel('Reflectance')
xlim([0.4 0.68])
ylim([0 1])
legend1 = sprintf('L = %0.2f nm', 1000*esti_L_lsqr);
legend2 = sprintf('L = %0.2f nm', 1000*theoric_L);
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

%%

%CODE TO SHOW AUTOMATICALLY CROPPED OXIDE AND SILICON IMAGES USING CUTOFF VALUE

figure(5)
subplot(1,2,1)
imshow(R_Tiff_Xnm_Crop) 
title('Red Oxide Image after cropping using cutoff value')
subplot(1,2,2)
imshow(R_Tiff_0nm_Crop)
title('Red Silicon Image after cropping using cutoff value')

figure(6)
subplot(1,2,1)
imshow(G_Tiff_Xnm_Crop) 
title('Green Oxide Image after cropping using cutoff value')
subplot(1,2,2)
imshow(G_Tiff_0nm_Crop)
title('Green Silicon Image after cropping using cutoff value')

figure(7)
subplot(1,2,1)
imshow(B_Tiff_Xnm_Crop)
title('Blue Oxide Image after cropping using cutoff value')
subplot(1,2,2)
imshow(B_Tiff_0nm_Crop)
title('Blue Silicon Image after cropping using cutoff value')

%CODE TO SHOW MANUALLY SELECTED OXIDE IMAGES

figure(8)
imshow(cropped_im_R_Xnm)
title('Cropped R Xnm ROI')
figure(9)
imshow(cropped_im_G_Xnm)
title('Cropped G Xnm ROI')
figure(10)
imshow(cropped_im_B_Xnm)
title('Cropped B Xnm ROI')

%%






