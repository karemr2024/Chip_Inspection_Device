clc; clear; clearvars;

% Import 'unknown' thickness images as RGB channels by prompting the user
[R_Tiff_Name,R_Tiff_Path] = uigetfile('*.tif','Red Image');
if isequal(R_Tiff_Name,0)
   disp('User selected Cancel');
else
   disp(['User selected ', fullfile(R_Tiff_Name,R_Tiff_Path)]);
   R_Tiff = imread(strcat(R_Tiff_Path, R_Tiff_Name));
end
[G_Tiff_Name,G_Tiff_Path] = uigetfile('*.tif','Green Image');
if isequal(G_Tiff_Name,0)
   disp('User selected Cancel');
else
   disp(['User selected ', fullfile(G_Tiff_Name,G_Tiff_Path)]);
   G_Tiff = imread(strcat(G_Tiff_Path, G_Tiff_Name));
end
[B_Tiff_Name,B_Tiff_Path] = uigetfile('*.tif','Blue Image');
if isequal(B_Tiff_Name,0)
   disp('User selected Cancel');
else
   disp(['User selected ', fullfile(B_Tiff_Name,B_Tiff_Path)]);
   B_Tiff = imread(strcat(B_Tiff_Path, B_Tiff_Name));
end
% Create RGB pixel masks according to camera sensor superpixel filter
% configuration 'rggb' and resolution 1080x1440.
R_Repeater = [1 0 ; 0 0];
G_Repeater = [0 1 ; 1 0];
B_Repeater = [0 0 ; 0 1];

R_Mask = uint16(repmat(R_Repeater,540,720));
G_Mask = uint16(repmat(G_Repeater,540,720));
B_Mask = uint16(repmat(B_Repeater,540,720));
% Multiply and add The original image under red light and the red mask 
% so that the R_Chan matrix represents only the light that
% red pixels of the sensor receive. Repeat for green and blue.
R_Chan = R_Tiff.*R_Mask;
G_Chan = G_Tiff.*G_Mask;
B_Chan = B_Tiff.*B_Mask;
% Combine channels to aquire new image.
RGB_Combo = R_Chan + G_Chan + B_Chan;
% Grayscale and demosaic the new image
Demo_Img = demosaic(RGB_Combo,'rggb');
figure(1)
f1 = imshow(Demo_Img);
imwrite(Demo_Img,'RGB_Combo_Demo.tif')

% Pixel leakage test:
% Import images of uniform mirror surface under different lighting conditions
% User selects image under no lighting
[DARK_Tiff_Name,DARK_Tiff_Path] = uigetfile('*.tif','Dark Image');
if isequal(DARK_Tiff_Name,0)
   disp('User selected Cancel');
else
   disp(['User selected ', fullfile(DARK_Tiff_Name,DARK_Tiff_Path)]);
   DARK_Tiff = imread(strcat(DARK_Tiff_Path, DARK_Tiff_Name));
end
% User selects image under red light
[RLight_Tiff_Name,RLight_Tiff_Path] = uigetfile('*.tif','Red Leak-Test Image');
if isequal(RLight_Tiff_Name,0)
   disp('User selected Cancel');
else
   disp(['User selected ', fullfile(RLight_Tiff_Name,RLight_Tiff_Path)]);
   RLight_Tiff = imread(strcat(RLight_Tiff_Path, RLight_Tiff_Name));
end
% User selects image under green light
[GLight_Tiff_Name,GLight_Tiff_Path] = uigetfile('*.tif','Green Leak-Test Image');
if isequal(GLight_Tiff_Name,0)
   disp('User selected Cancel');
else
   disp(['User selected ', fullfile(GLight_Tiff_Name,GLight_Tiff_Path)]);
   GLight_Tiff = imread(strcat(GLight_Tiff_Path, GLight_Tiff_Name));
end
% User selects image under blue light
[BLight_Tiff_Name,BLight_Tiff_Path] = uigetfile('*.tif','Blue Leak-Test Image');
if isequal(BLight_Tiff_Name,0)
   disp('User selected Cancel');
else
   disp(['User selected ', fullfile(BLight_Tiff_Name,BLight_Tiff_Path)]);
   BLight_Tiff = imread(strcat(BLight_Tiff_Path, BLight_Tiff_Name));
end
%%
% Lighting uniformity:
RLight_dem = demosaic(RLight_Tiff,'rggb');
GLight_dem = demosaic(GLight_Tiff,'rggb');
BLight_dem = demosaic(BLight_Tiff,'rggb');
figure(2)
fR = imshow(RLight_dem);
imwrite(rgb2gray(RLight_dem),'R_Light_dem.tif')
figure(3)
fG = imshow(GLight_dem);
imwrite(rgb2gray(GLight_dem),'G_Light_dem.tif')
figure(4)
fB = imshow(BLight_dem);
imwrite(rgb2gray(BLight_dem),'B_Light_dem.tif')

% Leakage: 
% Leakage of red light into blue bayer pixels
BLeak_for_RLight = (RLight_Tiff - DARK_Tiff).*B_Mask;
%Leakage of red light into green bayer pixels
GLeak_for_RLight = (RLight_Tiff- DARK_Tiff).*G_Mask;
% Leakage of green light into blue bayer pixels
BLeak_for_GLight = (GLight_Tiff - DARK_Tiff).*B_Mask;
% Leakage of green light into red bayer pixels
RLeak_for_GLight = (GLight_Tiff - DARK_Tiff).*R_Mask;
% Leakage of blue light into green bayer pixels
GLeak_for_BLight = (BLight_Tiff - DARK_Tiff).*G_Mask;
% Leakage of blue light into red bayer pixels
RLeak_for_BLight = (BLight_Tiff - DARK_Tiff).*R_Mask;
%%
R_im = RLight_Tiff - (RLeak_for_BLight + RLeak_for_GLight);
figure(5)
imshow(R_im)

G_im = GLight_Tiff - (GLeak_for_RLight + GLeak_for_BLight);
figure(6)
imshow(G_im)

B_im = BLight_Tiff - (BLeak_for_RLight + BLeak_for_GLight);
figure(7)
imshow(B_im)
%%
% Combination after leakage correction

RGB_Combo_Unleaked = (R_im.*R_Mask) + (G_im.*G_Mask) + (B_im.*B_Mask);

RGB_Combo_Unleaked_Demosaic = demosaic(RGB_Combo_Unleaked,'rggb');
figure(8)
imshow(RGB_Combo_Unleaked_Demosaic)
imwrite(RGB_Combo_Unleaked_Demosaic,'RGB_Combo_Unleaked_Demosaic.tif')

% Calculating thicknesses using reflectences from eren's code. 
% This part of the code has two major problems:
% 1. The way RefRed_at_XnmO2 (reflectance at 'unknown' SiO2 thickness) is
% calculated is prone to errors when G_im (image under green light,
% adjusted for leakage of other light at dark) has very small values
% compared to G_Chan (green channel of image of 'unknown' SiO2 thickness
% under green light).
% 2. The way mean reflectance values are calculated does not account for zero values that may
% occur in the pixels of interest.
%%
I_inR = (R_im.*R_Mask)./ .4028;
I_inG = (G_im.*G_Mask)./ .3720;
I_inB = (B_im.*B_Mask)./ .3448;

% I_inR = (255*((R_im.*R_Mask)./max(max(R_im.*R_Mask))))./ 40.28;
% I_inG = (255*((G_im.*G_Mask)./max(max((G_im.*G_Mask))))./ 37.20;
% I_inB = (255*((B_im.*G_Mask)./max(max(B_im.*G_Mask))))./ 34.48;

% Normalise intensity (Added on 18 MAY 2022)
I_inRn = 255*(I_inR./max(max(I_inR)));
I_inGn = 255*(I_inG./max(max(I_inG)));
I_inBn = 255*(I_inB./max(max(I_inB)));
%% 
RefRed_at_XnmO2 = (255*(R_Chan./max(max(R_Chan))))./I_inRn;
RefGre_at_XnmO2 = (255*(G_Chan./max(max(G_Chan))))./I_inGn;
RefBlu_at_XnmO2 = (255*(B_Chan./max(max(B_Chan))))./I_inBn;
%%
MeanRefRed_at_XnmO2 = mean(mean(nonzeros(RefRed_at_XnmO2)));
MeanRefGre_at_XnmO2 = mean(mean(nonzeros(RefGre_at_XnmO2)));
MeanRefBlu_at_XnmO2 = mean(mean(nonzeros(RefBlu_at_XnmO2)));

Ref_Fitter = [MeanRefRed_at_XnmO2 MeanRefGre_at_XnmO2 MeanRefBlu_at_XnmO2];
%%
% Leaktest Stats:
% Means of each nonzero pixel are taken for image taken under X colour (XLight_Tiff),
% means of X_im (image under X coloured light, adjusted for leakage of the other two colours of light at dark)
% are taken. For each colour:
% [(mean of non zero indeces of XLight_Tiff)-(mean of non zero indeces of X_im)]/(mean of non zero indeces of XLight_Tiff)
% is used to calculate leakage. This value is multiplied by 100 to get a
% percentage. 

R_Mean = mean(RLight_Tiff,2);
%R_Std = std(RLight_Tiff);
percentleak_R = 100*((mean(mean(nonzeros(RLight_Tiff)))-mean(mean(nonzeros(R_im))))/mean(mean(nonzeros(RLight_Tiff))));

G_Mean = mean(GLight_Tiff);
%G_Std = std(GLight_Tiff);
percentleak_G = 100*((mean(mean(nonzeros(GLight_Tiff)))-mean(mean(nonzeros(G_im))))/mean(mean(nonzeros(GLight_Tiff))));

B_Mean = mean(BLight_Tiff);
%B_Std = std(BLight_Tiff);
percentleak_B = 100*((mean(mean(nonzeros(BLight_Tiff)))-mean(mean(nonzeros(B_im))))/mean(mean(nonzeros(BLight_Tiff))));

fprintf('Leakege of green and blue light into red pixel is %f percent. \n',percentleak_R);
fprintf('Leakege of red and blue light into green pixel is %f percent. \n',percentleak_G);
fprintf('Leakege of red and green light into blue pixel is %f percent. \n',percentleak_B);

% Ellipsometry:
% https://www.jawoollam.com/resources/ellipsometry-tutorial/interaction-of-light-and-materials
% http://homes.nano.aau.dk/kp/Ellipsometry/main.pdf
% Thickness from reflectence
% https://opg.optica.org/ao/fulltext.cfm?uri=ao-48-5-985
% ehe