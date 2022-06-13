clc; close all; clearvars;

%Below are prompts for user to enter unknown SiO2 nm (120 nm for now)
%and 0 nm (Si) chip images.

[R_Tiff_Name,R_Tiff_Path] = uigetfile('*.tif','Red Image'); %Import Unknown Thickness: Red Image
if isequal(R_Tiff_Name,0)
   disp('User selected Cancel');
else
   disp(['User selected ', fullfile(R_Tiff_Name,R_Tiff_Path)]);
   R_Tiff_Xnm = imread(strcat(R_Tiff_Path, R_Tiff_Name));
end

[G_Tiff_Name,G_Tiff_Path] = uigetfile('*.tif','Green Image'); %Import Unknown Thickness: Green Image
if isequal(G_Tiff_Name,0)
   disp('User selected Cancel');
else
   disp(['User selected ', fullfile(G_Tiff_Name,G_Tiff_Path)]);
   G_Tiff_Xnm = imread(strcat(G_Tiff_Path, G_Tiff_Name));
end

[B_Tiff_Name,B_Tiff_Path] = uigetfile('*.tif','Blue Image'); %Import Unknown Thickness: Blue Image 
if isequal(B_Tiff_Name,0)
   disp('User selected Cancel');
else
   disp(['User selected ', fullfile(B_Tiff_Name,B_Tiff_Path)]);
   B_Tiff_Xnm = imread(strcat(B_Tiff_Path, B_Tiff_Name));
end

[DARK_Tiff_Name,DARK_Tiff_Path] = uigetfile('*.tif','Dark Image'); %Import 0 nm: Dark Image
if isequal(DARK_Tiff_Name,0)
   disp('User selected Cancel');
else
   disp(['User selected ', fullfile(DARK_Tiff_Name,DARK_Tiff_Path)]);
   Dark_Tiff_0nm = imread(strcat(DARK_Tiff_Path, DARK_Tiff_Name));
end

[RLight_Tiff_Name,RLight_Tiff_Path] = uigetfile('*.tif','Red Leak-Test Image'); %Import 0 nm: Red Image
if isequal(RLight_Tiff_Name,0)
   disp('User selected Cancel');
else
   disp(['User selected ', fullfile(RLight_Tiff_Name,RLight_Tiff_Path)]);
   R_Tiff_0nm = imread(strcat(RLight_Tiff_Path, RLight_Tiff_Name));
end

[GLight_Tiff_Name,GLight_Tiff_Path] = uigetfile('*.tif','Green Leak-Test Image'); %Import 0 nm: Green Image
if isequal(GLight_Tiff_Name,0)
   disp('User selected Cancel');
else
   disp(['User selected ', fullfile(GLight_Tiff_Name,GLight_Tiff_Path)]);
   G_Tiff_0nm = imread(strcat(GLight_Tiff_Path, GLight_Tiff_Name));
end

[BLight_Tiff_Name,BLight_Tiff_Path] = uigetfile('*.tif','Blue Leak-Test Image'); %Import 0 nm: Blue Image
if isequal(BLight_Tiff_Name,0)
   disp('User selected Cancel');
else
   disp(['User selected ', fullfile(BLight_Tiff_Name,BLight_Tiff_Path)]);
   B_Tiff_0nm = imread(strcat(BLight_Tiff_Path, BLight_Tiff_Name));
end

clearvars -except B_Tiff_0nm B_Tiff_Xnm Dark_Tiff_0nm G_Tiff_0nm G_Tiff_Xnm B_Tiff_0nm B_Tiff_Xnm R_Tiff_0nm R_Tiff_Xnm
%%

theoric_L = (input('Enter the thickness of your chip in nm (if unknown, enter UNKNOWN) : \n'))/1000;

load("Imaging_Data.mat")

R_Mask = uint16(repmat([1 0 ; 0 0],540,720)); % Create RGB pixel masks according 
G_Mask = uint16(repmat([0 1 ; 1 0],540,720)); % to camera sensor superpixel filter 
B_Mask = uint16(repmat([0 0 ; 0 1],540,720)); % configuration 'rggb' and resolution 1080x1440.

%For 0 nm image

R_Chan_0nm = demosaic(R_Tiff_0nm.*R_Mask,'rggb'); %Multiply Each R,G,B image with their individual
G_Chan_0nm = demosaic(G_Tiff_0nm.*G_Mask,'rggb'); %masks to take only the corresponding color light
B_Chan_0nm = demosaic(B_Tiff_0nm.*B_Mask,'rggb'); %on that pixel and create R,G,B channels

RGB_Combo_0nm = R_Chan_0nm + G_Chan_0nm + B_Chan_0nm; %Add R,G,B Channels to create a RGB combo image

%For Unknown thickness image

R_Chan_Xnm = demosaic(R_Tiff_Xnm.*R_Mask,'rggb'); %Multiply Each R,G,B image with their individual
G_Chan_Xnm = demosaic(G_Tiff_Xnm.*G_Mask,'rggb'); %masks to take only the corresponding color light
B_Chan_Xnm = demosaic(B_Tiff_Xnm.*B_Mask,'rggb'); %on that pixel and create R,G,B channels

RGB_Combo_Xnm = R_Chan_Xnm + G_Chan_Xnm + B_Chan_Xnm; %Add R,G,B Channels to create a RGB combo image

figure(1)
imagesc(R_Chan_0nm)
title('R Channel of 0 nm image')

figure(2)
imagesc(G_Chan_0nm)
title('G Channel of 0 nm image')

figure(3)
imagesc(B_Chan_0nm)
title('B Channel of 0 nm image')

figure(4)
imagesc(RGB_Combo_0nm)
title('RGB Image of 0 nm thickness')

figure(5)
imagesc(R_Chan_Xnm)
title('R Channel of unknown thickness image')

figure(6)
imagesc(G_Chan_Xnm)
title('G Channel of unknown thickness image')

figure(7)
imagesc(B_Chan_Xnm)
title('B Channel of unknown thickness image')

figure(8)
imagesc(RGB_Combo_Xnm)
title('RGB Image of unknown thickness')

%Leakage calculations for 0 nm image

RedLeak2GreenPix0nm = (R_Tiff_0nm - Dark_Tiff_0nm).*G_Mask;    %Red Leakage into Green Pixels
RedLeak2BluePix0nm = (R_Tiff_0nm - Dark_Tiff_0nm).*B_Mask;     %Red Leakage into Blue Pixels 
GreenLeak2RedPix0nm = (G_Tiff_0nm - Dark_Tiff_0nm).*R_Mask;    %Green Leakage into Red Pixels
GreenLeak2BluePix0nm = (G_Tiff_0nm - Dark_Tiff_0nm).*B_Mask;   %Green Leakage into Blue Pixels
BlueLeak2RedPix0nm = (B_Tiff_0nm - Dark_Tiff_0nm).*R_Mask;     %Blue Leakage into Red Pixels
BlueLeak2GreenPix0nm = (B_Tiff_0nm - Dark_Tiff_0nm).*G_Mask;   %Blue Leakage into Green Pixels

R_Chan_0nm_LeakCorrected = (R_Tiff_0nm - (BlueLeak2RedPix0nm + GreenLeak2RedPix0nm));
G_Chan_0nm_LeakCorrected = (G_Tiff_0nm - (RedLeak2GreenPix0nm + BlueLeak2GreenPix0nm));
B_Chan_0nm_LeakCorrected = (B_Tiff_0nm - (RedLeak2BluePix0nm + GreenLeak2BluePix0nm));

R_Chan_0nm_LeakCorrected_demo = demosaic(R_Chan_0nm_LeakCorrected,'rggb');
G_Chan_0nm_LeakCorrected_demo = demosaic(G_Chan_0nm_LeakCorrected,'rggb');
B_Chan_0nm_LeakCorrected_demo = demosaic(B_Chan_0nm_LeakCorrected,'rggb');

RGB_0nm_LeakCorrected = R_Chan_0nm_LeakCorrected_demo + G_Chan_0nm_LeakCorrected_demo + B_Chan_0nm_LeakCorrected_demo;

percentleak2R0nm = 100*((mean(mean(nonzeros(R_Tiff_0nm)))-mean(mean(nonzeros(R_Chan_0nm_LeakCorrected))))/mean(mean(nonzeros(R_Tiff_0nm))));
percentleak2G0nm = 100*((mean(mean(nonzeros(G_Tiff_0nm)))-mean(mean(nonzeros(G_Chan_0nm_LeakCorrected))))/mean(mean(nonzeros(G_Tiff_0nm))));
percentleak2B0nm = 100*((mean(mean(nonzeros(B_Tiff_0nm)))-mean(mean(nonzeros(B_Chan_0nm_LeakCorrected))))/mean(mean(nonzeros(B_Tiff_0nm))));

fprintf('Leakage of green and blue light into red pixel in 0 nm image is %f percent. \n',percentleak2R0nm);
fprintf('Leakage of red and blue light into green pixel in 0 nm image is %f percent. \n',percentleak2G0nm);
fprintf('Leakage of red and green light into blue pixel in 0 nm image is %f percent. \n',percentleak2B0nm);

% %Leakage calculations for unknown thickness image (DON'T HAVE 120 NM DARK IMAGE YET!!!)
% 
% RedLeak2GreenPixXnm = (R_Tiff_Xnm - Dark_Tiff_Xnm).*G_Mask;    %Red Leakage into Green Pixels
% RedLeak2BluePixXnm = (R_Tiff_Xnm - Dark_Tiff_Xnm).*B_Mask;     %Red Leakage into Blue Pixels 
% GreenLeak2RedPixXnm = (G_Tiff_Xnm - Dark_Tiff_Xnm).*R_Mask;    %Green Leakage into Red Pixels
% GreenLeak2BluePixXnm = (G_Tiff_Xnm - Dark_Tiff_Xnm).*B_Mask;   %Green Leakage into Blue Pixels
% BlueLeak2RedPixXnm = (B_Tiff_Xnm - Dark_Tiff_Xnm).*R_Mask;     %Blue Leakage into Red Pixels
% BlueLeak2GreenPixXnm = (B_Tiff_Xnm - Dark_Tiff_Xnm).*G_Mask;   %Blue Leakage into Green Pixels
% 
% R_Chan_Xnm_LeakCorrected = R_Tiff_Xnm - (BlueLeak2RedPixXnm + GreenLeak2RedPixXnm);
% G_Chan_Xnm_LeakCorrected = G_Tiff_Xnm - (RedLeak2GreenPixXnm + BlueLeak2GreenPixXnm);
% B_Chan_Xnm_LeakCorrected = B_Tiff_Xnm - (RedLeak2BluePixXnm + GreenLeak2BluePixXnm);
% 
% R_Chan_Xnm_LeakCorrected_demo = demosaic(R_Chan_Xnm_LeakCorrected,'rggb');
% G_Chan_Xnm_LeakCorrected_demo = demosaic(G_Chan_Xnm_LeakCorrected,'rggb');
% B_Chan_Xnm_LeakCorrected_demo = demosaic(B_Chan_Xnm_LeakCorrected,'rggb');
% 
% RGB_Xnm_LeakCorrected = R_Chan_Xnm_LeakCorrected_demo + G_Chan_Xnm_LeakCorrected_demo + B_Chan_Xnm_LeakCorrected_demo;
% 
% percentleak2RXnm = 100.*((mean(mean(R_Tiff_Xnm))-mean(mean(R_Chan_Xnm_LeakCorrected)))./mean(mean(R_Tiff_Xnm)));
% percentleak2GXnm = 100.*((mean(mean(G_Tiff_Xnm))-mean(mean(G_Chan_Xnm_LeakCorrected)))./mean(mean(G_Tiff_Xnm)));
% percentleak2BXnm = 100.*((mean(mean(B_Tiff_Xnm))-mean(mean(B_Chan_Xnm_LeakCorrected)))./mean(mean(B_Tiff_Xnm)));
% 
% fprintf('Leakage of green and blue light into red pixel in unknown thickness image is %f percent. \n',percentleak2RXnm);
% fprintf('Leakage of red and blue light into green pixel in unknown thickness image is %f percent. \n',percentleak2GXnm);
% fprintf('Leakage of red and green light into blue pixel in unknown thickness image is %f percent. \n',percentleak2BXnm);

figure(9)
imagesc(R_Chan_0nm_LeakCorrected_demo)
title('Leak Corrected R Channel of 0 nm image')

figure(10)
imagesc(G_Chan_0nm_LeakCorrected_demo)
title('Leak Corrected G Channel of 0 nm image')

figure(11)
imagesc(B_Chan_0nm_LeakCorrected_demo)
title('Leak Corrected B Channel of 0 nm image')

figure(12)
imagesc(RGB_0nm_LeakCorrected)
title('Leak Corrected RGB image of 0 nm thickness')

% figure(13)
% imagesc(R_Chan_Xnm_LeakCorrected_demo)
% title('Leak Corrected R Channel of 0 nm image')
% 
% figure(14)
% imagesc(G_Chan_0nm_LeakCorrected_demo)
% title('Leak Corrected G Channel of 0 nm image')
% 
% figure(15)
% imagesc(B_Chan_0nm_LeakCorrected_demo)
% title('Leak Corrected B Channel of 0 nm image')
% 
% figure(16)
% imagesc(RGB_0nm_LeakCorrected)
% title('Leak Corrected RGB image of 0 nm thickness')

%Combination for 0 nm after leakage correction

RGB_Combo_Unleaked_0nm = (R_Chan_0nm_LeakCorrected.*R_Mask) + (G_Chan_0nm_LeakCorrected.*G_Mask) + (B_Chan_0nm_LeakCorrected.*B_Mask);
RGB_Combo_Unleaked_0nm_Demosaic = demosaic(RGB_Combo_Unleaked_0nm,'rggb');

%Reflectances at 0 and 120 nm

Ref_at_0nm_R = Gamma(valtoindex_lambda(cw_r),valtoindex_L(0));
Ref_at_0nm_G = Gamma(valtoindex_lambda(cw_g),valtoindex_L(0));
Ref_at_0nm_B = Gamma(valtoindex_lambda(cw_b),valtoindex_L(0));
Ref_at_120nm_R = Gamma(valtoindex_lambda(cw_r),valtoindex_L(0.12));
Ref_at_120nm_G = Gamma(valtoindex_lambda(cw_g),valtoindex_L(0.12));
Ref_at_120nm_B = Gamma(valtoindex_lambda(cw_b),valtoindex_L(0.12));


[R_Tiff_0nm_Crop,x_crop_0nm,y_crop_0nm] = AreaSelection(R_Tiff_0nm);
G_Tiff_0nm_Crop = G_Tiff_0nm(x_crop_0nm:x_crop_0nm+63,y_crop_0nm:y_crop_0nm+63);
B_Tiff_0nm_Crop = B_Tiff_0nm(x_crop_0nm:x_crop_0nm+63,y_crop_0nm:y_crop_0nm+63);

[R_Tiff_Xnm_Crop,x_crop_Xnm,y_crop_Xnm] = AreaSelection(R_Tiff_Xnm);
G_Tiff_Xnm_Crop = G_Tiff_Xnm(x_crop_Xnm:x_crop_Xnm+63,y_crop_Xnm:y_crop_Xnm+63);
B_Tiff_Xnm_Crop = B_Tiff_Xnm(x_crop_Xnm:x_crop_Xnm+63,y_crop_Xnm:y_crop_Xnm+63);

R_Mask_Xnm_Crop = double(R_Tiff_Xnm_Crop > 8000);

if  G_Tiff_Xnm_Crop(1,2) < G_Tiff_Xnm_Crop(1,1) && G_Tiff_Xnm_Crop(1,2) < G_Tiff_Xnm_Crop(2,1) && G_Tiff_Xnm_Crop(1,2) < G_Tiff_Xnm_Crop(2,2) 
    G_Mask_Xnm_Crop = double(repmat([1 0; 0 1],32,32));
end

B_Mask_Xnm_Crop = double(B_Tiff_Xnm_Crop > 10000);

%Intensity In calculation (calibration) I_in = I_out_Si / Reflectance

I_inR = round(double(R_Tiff_0nm_Crop) ./ Ref_at_0nm_R);
I_inG = round(double(G_Tiff_0nm_Crop) ./ Ref_at_0nm_G);
I_inB = round(double(B_Tiff_0nm_Crop) ./ Ref_at_0nm_B);

%Reflectance from oxide calculation: Reflectance_Oxide = I_Out_SiO2 / I_in
%Multiplied with masks to eliminate other pixels and thus divide by zeros.

RefRed_at_XnmO2 = (double(R_Tiff_Xnm_Crop)./I_inR).*double(R_Mask_Xnm_Crop);
RefGre_at_XnmO2 = (double(G_Tiff_Xnm_Crop)./I_inG).*double(G_Mask_Xnm_Crop);
RefBlu_at_XnmO2 = (double(B_Tiff_Xnm_Crop)./I_inB).*double(B_Mask_Xnm_Crop);

%Take the average of Reflectance values in R,G,B only pixels

MeanRefRed_at_XnmO2 = mean(mean(nonzeros(RefRed_at_XnmO2)));
MeanRefGre_at_XnmO2 = mean(mean(nonzeros(RefGre_at_XnmO2)));
MeanRefBlu_at_XnmO2 = mean(mean(nonzeros(RefBlu_at_XnmO2)));

%use reftocurve function to estimate thickness from these reflectance values
esti_L = reftocurve(MeanRefRed_at_XnmO2,MeanRefGre_at_XnmO2,MeanRefBlu_at_XnmO2);
esti_L_lsqr = reftocurve_lsqr(MeanRefRed_at_XnmO2,MeanRefGre_at_XnmO2,MeanRefBlu_at_XnmO2);

if isnumeric(theoric_L)
percent_error = 100*((abs(esti_L_lsqr-theoric_L))/theoric_L);
end

figure(13)

hold on
plot(lambda,Gamma(:,valtoindex_L(abs(esti_L_lsqr))),'m--','LineWidth',2)           %Estimated nm thickness
plot(lambda,Gamma(:,valtoindex_L(0.12)),'c--','LineWidth',2)                       %120 nm thickness

plot(cw_r,MeanRefRed_at_XnmO2,'r.','MarkerSize',30)
plot(cw_g,MeanRefGre_at_XnmO2,'g.','MarkerSize',30)
plot(cw_b,MeanRefBlu_at_XnmO2,'b.','MarkerSize',30)

xlabel('lambda (\mum)')
ylabel('Reflectance')
xlim([0.4 0.68])
ylim([0 1])
legend1 = sprintf('L = %0.2f nm', 1000*esti_L_lsqr);
legend2 = sprintf('L = 120 nm');
legend(legend1, legend2,'Ref at R','Ref at G','Ref at B','location','bestoutside')
fprintf('Actual thickness is %f nm \n',theoric_L*1000)
fprintf('Estimated thickness is %f nm \n',esti_L_lsqr*1000)
fprintf('Percent Error is %f%',percent_error)






