clc; close all; clearvars -except Gamma_Si;

%Below are prompts for user to enter unknown SiO2 nm (120 nm for now)
%and 0 nm (Si) chip images.

%UNCOMMENT STARTING HERE

% [R_Tiff_Name,R_Tiff_Path] = uigetfile('*.tif','Red Image'); %Import Unknown Thickness: Red Image
% if isequal(R_Tiff_Name,0)
%    disp('User selected Cancel');
% else
%    disp(['User selected ', fullfile(R_Tiff_Name,R_Tiff_Path)]);
%    R_Tiff = imread(strcat(R_Tiff_Path, R_Tiff_Name));
% end
% 
% [G_Tiff_Name,G_Tiff_Path] = uigetfile('*.tif','Green Image'); %Import Unknown Thickness: Green Image
% if isequal(G_Tiff_Name,0)
%    disp('User selected Cancel');
% else
%    disp(['User selected ', fullfile(G_Tiff_Name,G_Tiff_Path)]);
%    G_Tiff = imread(strcat(G_Tiff_Path, G_Tiff_Name));
% end
% 
% [B_Tiff_Name,B_Tiff_Path] = uigetfile('*.tif','Blue Image'); %Import Unknown Thickness: Blue Image 
% if isequal(B_Tiff_Name,0)
%    disp('User selected Cancel');
% else
%    disp(['User selected ', fullfile(B_Tiff_Name,B_Tiff_Path)]);
%    B_Tiff = imread(strcat(B_Tiff_Path, B_Tiff_Name));
% end

%UNCOMMENT ENDS

tiff_info_R = imfinfo('StackR.tif'); % return tiff structure, one element per image
tiff_stack_R = imread('StackR.tif', 1) ; % read in first image
tiff_info_G = imfinfo('StackG.tif'); % return tiff structure, one element per image
tiff_stack_G = imread('StackG.tif', 1) ; % read in first image
tiff_info_B = imfinfo('StackG.tif'); % return tiff structure, one element per image
tiff_stack_B = imread('StackG.tif', 1) ; % read in first image
%concatenate each successive tiff to tiff_stack
for i = 2 : size(tiff_info_R, 1)
    temp_tiff_R = imread('StackR.tif', i);
    tiff_stack_R = cat(3 , tiff_stack_R, temp_tiff_R);
    temp_tiff_G = imread('StackG.tif', i);
    tiff_stack_G = cat(3 , tiff_stack_G, temp_tiff_G);
    temp_tiff_B = imread('StackB.tif', i);
    tiff_stack_B = cat(3 , tiff_stack_B, temp_tiff_B);
end

tiff_stack_R = double(tiff_stack_R);
tiff_stack_G = double(tiff_stack_G);
tiff_stack_B = double(tiff_stack_B);

firstim = 31;
lastim = 60;

sample_R = tiff_stack_R(:,:,firstim:lastim);
sample_G = tiff_stack_G(:,:,firstim:lastim);
sample_B = tiff_stack_B(:,:,firstim:lastim);
%
[rows,cols,~] = size(tiff_stack_R);

sampleavg_R = uint16(zeros(rows,cols));
sampleavg_G = uint16(zeros(rows,cols));
sampleavg_B = uint16(zeros(rows,cols));
%
sampleavg_R = sum(sample_R,3);
sampleavg_G = sum(sample_G,3);
sampleavg_B = sum(sample_B,3);
%
R_Tiff = uint16(sampleavg_R./(lastim-firstim+1));
G_Tiff = uint16(sampleavg_G./(lastim-firstim+1));
B_Tiff = uint16(sampleavg_B./(lastim-firstim+1));


%
clearvars -except B_Tiff G_Tiff R_Tiff 
%%
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

[cropped_im_R_0nm_in,x_crop_0nm_in,y_crop_0nm_in] = AreaSelection_Circle_SiIn(cropped_im_R_0nm_out,x_crop_0nm_out,y_crop_0nm_out);
cropped_im_G_0nm_in = AreaSelection_Circle_Mod_SiIn(cropped_im_G_0nm_out,x_crop_0nm_in,y_crop_0nm_in);
cropped_im_B_0nm_in = AreaSelection_Circle_Mod_SiIn(cropped_im_B_0nm_out,x_crop_0nm_in,y_crop_0nm_in);

[cropped_im_R_Xnm,x_crop_Xnm,y_crop_Xnm] = AreaSelection_Circle_Ox(cropped_im_R_0nm_in);
cropped_im_G_Xnm = AreaSelection_Circle_Mod(cropped_im_G_0nm_in,x_crop_Xnm,y_crop_Xnm);
cropped_im_B_Xnm = AreaSelection_Circle_Mod(cropped_im_B_0nm_in,x_crop_Xnm,y_crop_Xnm);


thin_ring_R = uint16(cropped_im_R_0nm_out) - uint16(cropped_im_R_0nm_in);
thin_ring_G = uint16(cropped_im_G_0nm_out) - uint16(cropped_im_G_0nm_in);
thin_ring_B = uint16(cropped_im_B_0nm_out) - uint16(cropped_im_B_0nm_in);


%%

nonzeros_R_0nm = nonzeros(thin_ring_R);
nonzeros_G_0nm = nonzeros(thin_ring_G);
nonzeros_B_0nm = nonzeros(thin_ring_B);

nonzeros_R_Xnm = nonzeros(cropped_im_R_Xnm);
nonzeros_G_Xnm = nonzeros(cropped_im_G_Xnm);
nonzeros_B_Xnm = nonzeros(cropped_im_B_Xnm);
%%
if numel(nonzeros_R_0nm) > numel(nonzeros_R_Xnm)
    nonzeros_R_0nm = nonzeros_R_0nm(1:numel(nonzeros_R_Xnm));
elseif numel(nonzeros_R_0nm) < numel(nonzeros_R_Xnm)
    nonzeros_R_Xnm = nonzeros_R_Xnm(1:numel(nonzeros_R_0nm)); 
end

if numel(nonzeros_G_0nm) > numel(nonzeros_G_Xnm)
    nonzeros_G_0nm = nonzeros_G_0nm(1:numel(nonzeros_G_Xnm));
elseif numel(nonzeros_G_0nm) < numel(nonzeros_G_Xnm)
    nonzeros_G_Xnm = nonzeros_G_Xnm(1:numel(nonzeros_G_0nm)); 
end

if numel(nonzeros_B_0nm) > numel(nonzeros_B_Xnm)
    nonzeros_B_0nm = nonzeros_B_0nm(1:numel(nonzeros_B_Xnm));
elseif numel(nonzeros_B_0nm) < numel(nonzeros_B_Xnm)
    nonzeros_B_Xnm = nonzeros_B_Xnm(1:numel(nonzeros_B_0nm)); 
end


%Intensity In calculation (calibration) I_in = I_out_Si / Reflectance

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

load("Gamma_Si.mat")
load("Gamma_R.mat")
load("Gamma_G.mat")
load("Gamma_B.mat")
%%
figure(15)

hold on
plot(lambda,Gamma_B(:,valtoindex_L(abs(esti_L_lsqr))),'r--','LineWidth',2)%Estimated nm thickness
plot(lambda,Gamma_G(:,valtoindex_L(abs(esti_L_lsqr))),'g--','LineWidth',2)
plot(lambda,Gamma_R(:,valtoindex_L(abs(esti_L_lsqr))),'b--','LineWidth',2)
% plot(lambda,Gamma(:,valtoindex_L(abs(esti_L_lsqr))),'b--','LineWidth',2)

plot(lambda,Gamma_B(:,valtoindex_L(theoric_L)),'c--','LineWidth',2)
plot(lambda,Gamma_G(:,valtoindex_L(theoric_L)),'c--','LineWidth',2)
plot(lambda,Gamma_R(:,valtoindex_L(theoric_L)),'c--','LineWidth',2)%120 nm thickness
% plot(lambda,Gamma(:,valtoindex_L(theoric_L)),'c--','LineWidth',2)
% 
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
title("Actual Thickness is " + theoric_L*1000 + " nm, Estimated Thickness is " + esti_L_lsqr*1000 + " nm")

if isnumeric(theoric_L)
percent_error = 100*((abs(esti_L_lsqr-theoric_L))/theoric_L);
fprintf('Actual thickness is %f nm \n',theoric_L*1000)
fprintf('Estimated thickness is %f nm \n',esti_L_lsqr*1000)
fprintf("Percent error is %f", percent_error)
elseif theoric_L == 'UNKNOWN'
fprintf('Estimated thickness is %f nm \n',esti_L_lsqr*1000)
disp('Percent Error is not known')
end

%Thickness mapping



%%SUPPLEMENTARY IMAGES AND FIGURES

%%
%%CODE TO SHOW MANUALLY SELECTED OXIDE IMAGES
%%
% figure(8)
% imshow(cropped_im_R_Xnm)
% title('Cropped R Xnm ROI')
% figure(9)
% imshow(cropped_im_G_Xnm)
% title('Cropped G Xnm ROI')
% figure(10)
% imshow(cropped_im_B_Xnm)
% title('Cropped B Xnm ROI')
%%
%CODE TO SHOW R G B RINGS
%%
figure(31)
imagesc(thin_ring_R)
colormap
colorbar
caxis([min(min(nonzeros(thin_ring_R))) max(max(nonzeros(thin_ring_R)))]);
title('R Ring ROI')
figure(32)
imagesc(thin_ring_G)
colormap
colorbar
caxis([min(min(nonzeros(thin_ring_G))) max(max(nonzeros(thin_ring_G)))]);
title('G Ring ROI')
figure(33)
imagesc(thin_ring_B)
colormap
colorbar
caxis([double(min(min(nonzeros(thin_ring_B)))) double(max(max(nonzeros(thin_ring_B))))]);
title('B Ring ROI')
ans = 1;    
%%






