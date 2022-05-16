clc; clear all; clearvars;

% Import images as RGB channels
R_Tiff = imread("Enter R Image Name");
G_Tiff = imread("Enter G Image Name");
B_Tiff = imread("Enter B Image Name");
% Create RGB pixel masks
R_Repeater = [1 0 ; 0 0];
G_Repeater = [0 1 ; 1 0];
B_Repeater = [0 0 ; 0 1];

R_Mask = uint16(repmat(R_Repeater,540,720));
G_Mask = uint16(repmat(G_Repeater,540,720));
B_Mask = uint16(repmat(B_Repeater,540,720));
% Multiply and add
R_Chan = R_Tiff.*R_Mask;
G_Chan = G_Tiff.*G_Mask;
B_Chan = B_Tiff.*B_Mask;

RGB_Combo = R_Chan + G_Chan + B_Chan;
% Grayscale and demosaic
Demo_Img = demosaic(RGB_Combo,'rggb');
figure(1)
f1 = imshow(Demo_Img);
imwrite(Demo_Img,'RGB_Combo_Demo.tif')

% Pixel leakage test:
% Import images of uniform mirror surface under different lighting conditions
DARK_Tiff = imread("Insert Dark Leak Test Image");
RLight_Tiff = imread("Insert Leak R Test Image");
GLight_Tiff = imread("Insert Leak G Test Image");
BLight_Tiff = imread("Insert Leak B Test Image");

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


% R_Mean = mean(RLight_Tiff);
% R_Std = std(RLight_Tiff);
% 
% G_Mean = mean(GLight_Tiff);
% G_Std = std(GLight_Tiff);
% 
% B_Mean = mean(BLight_Tiff);
% B_Std = std(BLight_Tiff);

% Leakage: 
% Leakage of red light into blue bayer pixels
BLeak_for_RLight = (RLight_Tiff - DARK_Tiff).*B_Mask;
%Leakage of red light into green bayer pixels
GLeak_for_RLight = (RLight_Tiff- DARK_Tiff).*G_Mask;
% Leakage of green light into blue bayer pixels
BLeak_for_GLight = (GLight_Tiff - DARK_Tiff).*B_Mask;
% Leakage of green light into red bayer pixels
RLeak_for_GLight = (GLight_Tiff- DARK_Tiff).*R_Mask;
% Leakage of blue light into green bayer pixels
GLeak_for_BLight = (BLight_Tiff - DARK_Tiff).*G_Mask;
% Leakage of blue light into red bayer pixels
RLeak_for_BLight = (BLight_Tiff - DARK_Tiff).*R_Mask;

% Image Correction (EREN)

% R_im = RLight_dem - RLeak_for_BLight -RLeak_for_GLight;
% figure(5)
% imshow(R_im)
% 
% G_im = GLight_dem - GLeak_for_RLight -GLeak_for_BLight;
% figure(6)
% imshow(G_im)
% 
% B_im = BLight_dem - BLeak_for_RLight -BLeak_for_GLight;
% figure(7)
% imshow(B_im)

R_im = RLight_Tiff - RLeak_for_BLight -RLeak_for_GLight;
figure(5)
imshow(R_im)

G_im = GLight_Tiff - GLeak_for_RLight -GLeak_for_BLight;
figure(6)
imshow(G_im)

B_im = BLight_Tiff - BLeak_for_RLight -BLeak_for_GLight;
figure(7)
imshow(B_im)

% Combination after leakage correction

RGB_Combo_Unleaked = R_im + G_im + B_im;

% I_inR = R_im./ RefMatR
% I_inG = G_im./ RefMatG
% I_inB = B_im./ RefMatB

RGB_Combo_Unleaked_Demosaic = demosaic(RGB_Combo_Unleaked,'rggb');
figure(8)
imshow(RGB_Combo_Unleaked_Demosaic)
imwrite(RGB_Combo_Unleaked_Demosaic,'RGB_Combo_Unleaked_Demosaic.tif')

% Ellipsometry:
% https://www.jawoollam.com/resources/ellipsometry-tutorial/interaction-of-light-and-materials
% http://homes.nano.aau.dk/kp/Ellipsometry/main.pdf
% Thickness from reflectence
% https://opg.optica.org/ao/fulltext.cfm?uri=ao-48-5-985
