%% START

tic
clc; clearvars; close all;
fprintf('Beginning to run %s.m ...\n', mfilename);
%%
addpath(genpath(pwd))

%% START NEW IMAGE LOAD

matB = dir('bim')
matG = dir('gim')
matR = dir('rim')

%% ENTER EXPERIMENT DETAILS AND CREATE FOLDER

% Chip_no = input('Chip #: ',"s")
% Exposure = input('Exposure time (ms): ',"s")
% FPS = input('FPS: ',"s")
% ImgCount = input('# of images: ' ,"s")
% Camera = 'BFS-U3-17S7M-C' 

%%

% foldername =  append('/experiments/results/CHIP',Chip_no,'_',Exposure,'ms_',FPS,'FPS_',ImgCount,'_IMG',Camera)
% projectdir = 'C:\Tepegoz\iRiS_Kinetics_Github\experiments\results';
% subfolder_name=char(append('/CHIP',Chip_no,'_',Exposure,'ms_',FPS,'FPS',ImgCount,'_IMG',Camera,'/CHIP_',Chip_no,'_',Exposure,'ms_',FPS,'FPS',ImgCount,'_IMG',Camera));
% mkdir(fullfile(projectdir, subfolder_name));

%%
matB = (matB(3:end));
matG = (matG(3:end));
matR = (matR(3:end));
%%
cellB = struct2cell(matB);
cellG = struct2cell(matG);
cellR = struct2cell(matR);
%%

imsB = cellB(1,:);
imsG = cellG(1,:);
imsR = cellR(1,:);

%%
[rows, cols] = size(cell2mat(imsR(1)));
nfiles = length(imsR);
%%

for i = 1:nfiles
B_Stack_pre(i) = load(imsB{1,i});
G_Stack_pre(i) = load(imsG{1,i});
R_Stack_pre(i) = load(imsR{1,i});
end

%% END NEW IMAGE LOAD

%% Below are prompts for user to input unknown SiO2 nm images in .mat format
% 
% [R_Tiff_Name,R_Tiff_Path] = uigetfile('*.mat','Red Stack'); %Import Unknown Thickness: Red Image
% if isequal(R_Tiff_Name,0)
%    disp('User selected Cancel');
% else
%    disp(['User selected ', fullfile(R_Tiff_Name,R_Tiff_Path)]);
%    R_Stack_pre = load(strcat(R_Tiff_Path, R_Tiff_Name));
% end
% 
% [G_Tiff_Name,G_Tiff_Path] = uigetfile('*.mat','Green Stack'); %Import Unknown Thickness: Green Image
% if isequal(G_Tiff_Name,0)
%    disp('User selected Cancel');
% else
%    disp(['User selected ', fullfile(G_Tiff_Name,G_Tiff_Path)]);
%    G_Stack_pre = load(strcat(G_Tiff_Path, G_Tiff_Name));
% end
% 
% [B_Tiff_Name,B_Tiff_Path] = uigetfile('*.mat','Blue Stack'); %Import Unknown Thickness: Blue Image
% if isequal(B_Tiff_Name,0)
%    disp('User selected Cancel');
% else
%    disp(['User selected ', fullfile(B_Tiff_Name,B_Tiff_Path)]);
%    B_Stack_pre = load(strcat(B_Tiff_Path, B_Tiff_Name));
% end

%% Prompt user about the thickness of the chip

load("Simulation_Data.mat")
% prompt_NorS = input('Do you know the thickness of your chip? (Y/N)\n', "s");
% if prompt_NorS == 'Y'
% theoric_L = (input('Enter the thickness of your chip in nm: \n'))/1000;
% else
% theoric_L = 'UNKNOWN';
% end

%% Convert RGB stacks from struct to cell

R_Stack = struct2cell(R_Stack_pre);
G_Stack = struct2cell(G_Stack_pre);
B_Stack = struct2cell(B_Stack_pre);

%% Separate images and timestamps

%% UN/COMMENT FOR OLD VERSION
% Rims = R_Stack(1:2:end);
% Rtimes =R_Stack(2:2:end);
% Gims = G_Stack(1:2:end);
% Gtimes = G_Stack(2:2:end);
% Bims = B_Stack(1:2:end);
% Btimes = B_Stack(2:2:end);
%% UN/COMMENT ENDS

%% UN/COMMENT FOR NEW VERSION
Rims = R_Stack;
Gims = G_Stack;
Bims = B_Stack;
%% UN/COMMENT ENDS

%% Get size

[rows, cols] = size(cell2mat(Rims(1)));
nfiles = length(Rims);

%% Convert RGB stacks to matrices

% tiff_stack_R = uint16(zeros(rows,cols,nfiles));
% tiff_stack_G = uint16(zeros(rows,cols,nfiles));
% tiff_stack_B = uint16(zeros(rows,cols,nfiles));
%%

for i = 1:nfiles
   tiff_stack_R(:,:,i) = cell2mat(Rims(i));
   tiff_stack_G(:,:,i) = cell2mat(Gims(i));
   tiff_stack_B(:,:,i) = cell2mat(Bims(i));
end
%%

%% Prompt user to select circular ROI and crop the same ROIs from each images

%% Display whole image and select ROI for spot

[R_Crop_Bothnm(:,:,1),x_crop_Bothnm,y_crop_Bothnm] = AreaSelection_Circle(tiff_stack_R(:,:,1));
G_Crop_Bothnm(:,:,1) = AreaSelection_Circle_Mod(tiff_stack_G(:,:,1),x_crop_Bothnm,y_crop_Bothnm);
B_Crop_Bothnm(:,:,1) = AreaSelection_Circle_Mod(tiff_stack_B(:,:,1),x_crop_Bothnm,y_crop_Bothnm);

for i = 2:nfiles
    R_Crop_Bothnm(:,:,i) = AreaSelection_Circle_Mod(tiff_stack_R(:,:,i),x_crop_Bothnm,y_crop_Bothnm);
    G_Crop_Bothnm(:,:,i) = AreaSelection_Circle_Mod(tiff_stack_G(:,:,i),x_crop_Bothnm,y_crop_Bothnm);
    B_Crop_Bothnm(:,:,i) = AreaSelection_Circle_Mod(tiff_stack_B(:,:,i),x_crop_Bothnm,y_crop_Bothnm);
end
%%

%% Display ROI selected spot image and select ROI for outer Silicon ring

[cropped_im_R_0nm_out(:,:,1),x_crop_0nm_out,y_crop_0nm_out,Center] = AreaSelection_Circle_SiOut(R_Crop_Bothnm(:,:,1));
cropped_im_G_0nm_out(:,:,1) = AreaSelection_Circle_Mod(G_Crop_Bothnm(:,:,1),x_crop_0nm_out,y_crop_0nm_out);
cropped_im_B_0nm_out(:,:,1) = AreaSelection_Circle_Mod(B_Crop_Bothnm(:,:,1),x_crop_0nm_out,y_crop_0nm_out);

for i = 2:nfiles
    cropped_im_R_0nm_out(:,:,i) = AreaSelection_Circle_Mod(R_Crop_Bothnm(:,:,i),x_crop_0nm_out,y_crop_0nm_out);
    cropped_im_G_0nm_out(:,:,i) = AreaSelection_Circle_Mod(G_Crop_Bothnm(:,:,i),x_crop_0nm_out,y_crop_0nm_out);
    cropped_im_B_0nm_out(:,:,i) = AreaSelection_Circle_Mod(B_Crop_Bothnm(:,:,i),x_crop_0nm_out,y_crop_0nm_out);
end

%% Display ROI selected outer ring image and select ROI for inner Silicon ring
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
[cropped_im_R_0nm_in(:,:,1),x_crop_0nm_in,y_crop_0nm_in] = AreaSelection_Circle_SiIn(R_Crop_Bothnm(:,:,1),x_crop_0nm_out,y_crop_0nm_out,Center);
cropped_im_G_0nm_in(:,:,1) = AreaSelection_Circle_Mod_SiIn(G_Crop_Bothnm(:,:,1),x_crop_0nm_in,y_crop_0nm_in);
cropped_im_B_0nm_in(:,:,1) = AreaSelection_Circle_Mod_SiIn(B_Crop_Bothnm(:,:,1),x_crop_0nm_in,y_crop_0nm_in);

for i = 2:nfiles
    cropped_im_R_0nm_in(:,:,i) = AreaSelection_Circle_Mod_SiIn(R_Crop_Bothnm(:,:,i),x_crop_0nm_in,y_crop_0nm_in);
    cropped_im_G_0nm_in(:,:,i) = AreaSelection_Circle_Mod_SiIn(G_Crop_Bothnm(:,:,i),x_crop_0nm_in,y_crop_0nm_in);
    cropped_im_B_0nm_in(:,:,i) = AreaSelection_Circle_Mod_SiIn(B_Crop_Bothnm(:,:,i),x_crop_0nm_in,y_crop_0nm_in);
end

%% Display ROI selected inner ring image and select ROI for Silicon Oxide

[cropped_im_R_Xnm(:,:,1),x_crop_Xnm,y_crop_Xnm] = AreaSelection_Circle_Ox(cropped_im_R_0nm_in(:,:,1));
cropped_im_G_Xnm(:,:,1) = AreaSelection_Circle_Mod(cropped_im_G_0nm_in(:,:,1),x_crop_Xnm,y_crop_Xnm);
cropped_im_B_Xnm(:,:,1) = AreaSelection_Circle_Mod(cropped_im_B_0nm_in(:,:,1),x_crop_Xnm,y_crop_Xnm);

for i = 2:nfiles
    cropped_im_R_Xnm(:,:,i) = AreaSelection_Circle_Mod(cropped_im_R_0nm_in(:,:,i),x_crop_Xnm,y_crop_Xnm);
    cropped_im_G_Xnm(:,:,i) = AreaSelection_Circle_Mod(cropped_im_G_0nm_in(:,:,i),x_crop_Xnm,y_crop_Xnm);
    cropped_im_B_Xnm(:,:,i) = AreaSelection_Circle_Mod(cropped_im_B_0nm_in(:,:,i),x_crop_Xnm,y_crop_Xnm);
end

%% Get the silicon rings by subtracting inner ring from outer ring for RGB

thin_ring_R = uint16(cropped_im_R_0nm_out) - uint16(cropped_im_R_0nm_in);
thin_ring_G = uint16(cropped_im_G_0nm_out) - uint16(cropped_im_G_0nm_in);
thin_ring_B = uint16(cropped_im_B_0nm_out) - uint16(cropped_im_B_0nm_in);

%%


%% Get nonzeros of one

nonzeros_R_0nm = nonzeros(thin_ring_R(:,:,1));
nonzeros_R_Xnm = nonzeros(cropped_im_R_Xnm(:,:,1));

%% Get size

[rowsnonzeros0nm, ~] = size(nonzeros_R_0nm);
[rowsnonzerosXnm, ~] = size(nonzeros_R_Xnm);

%% Preallocate

nonzeros_R_0nm = zeros(rowsnonzeros0nm, nfiles);
nonzeros_G_0nm = zeros(rowsnonzeros0nm, nfiles);
nonzeros_B_0nm = zeros(rowsnonzeros0nm, nfiles);

nonzeros_R_Xnm = zeros(rowsnonzerosXnm, nfiles);
nonzeros_G_Xnm = zeros(rowsnonzerosXnm, nfiles);
nonzeros_B_Xnm = zeros(rowsnonzerosXnm, nfiles);

%% Get nonzeros of each RGB Silicon ring

for i = 1:nfiles
nonzeros_R_0nm(:,i) = nonzeros(thin_ring_R(:,:,i));
nonzeros_G_0nm(:,i) = nonzeros(thin_ring_G(:,:,i));
nonzeros_B_0nm(:,i) = nonzeros(thin_ring_B(:,:,i));
end

%% Get nonzeros of each RGB Silicon Oxide ROI

for i = 1:nfiles
nonzeros_R_Xnm(:,i) = nonzeros(cropped_im_R_Xnm(:,:,i))';
nonzeros_G_Xnm(:,i) = nonzeros(cropped_im_G_Xnm(:,:,i))';
nonzeros_B_Xnm(:,i) = nonzeros(cropped_im_B_Xnm(:,:,i))';
end

%% Make the number of elements in silicon and oxide equal

% if numel(nonzeros_R_0nm) > numel(nonzeros_R_Xnm)
%     nonzeros_R_0nm = nonzeros_R_0nm(1:rowsnonzerosXnm,:);
% elseif numel(nonzeros_R_0nm) < numel(nonzeros_R_Xnm)
%     nonzeros_R_Xnm = nonzeros_R_Xnm(1:rowsnonzeros0nm,:); 
% end
% 
% if numel(nonzeros_G_0nm) > numel(nonzeros_G_Xnm)
%     nonzeros_G_0nm = nonzeros_G_0nm(1:rowsnonzerosXnm,:);
% elseif numel(nonzeros_G_0nm) < numel(nonzeros_G_Xnm)
%     nonzeros_G_Xnm = nonzeros_G_Xnm(1:rowsnonzeros0nm,:); 
% end
% 
% if numel(nonzeros_B_0nm) > numel(nonzeros_B_Xnm)
%     nonzeros_B_0nm = nonzeros_B_0nm(1:rowsnonzerosXnm,:);
% elseif numel(nonzeros_B_0nm) < numel(nonzeros_B_Xnm)
%     nonzeros_B_Xnm = nonzeros_B_Xnm(1:rowsnonzeros0nm,:); 
% end

%% INTENSITY v TIME GRAPHS AT SILICON ONLY (THIN RING)

%% Get the mean of the silicon values

meanR_0nm = mean(nonzeros_R_0nm);
meanG_0nm = mean(nonzeros_G_0nm);
meanB_0nm = mean(nonzeros_B_0nm);

%% Get the mean of the silicon oxide values

meanR_Xnm = mean(nonzeros_R_Xnm);
meanG_Xnm = mean(nonzeros_G_Xnm);
meanB_Xnm = mean(nonzeros_B_Xnm);

%% OPTIONAL FILTERING

%  for i = 1:nfiles
%      if avgR(i) > 1.05*meanavgR | avgR(i) < 0.95*meanavgR
%          avgR(i) = meanavgR;
%      end
%  end
% 
%  for i = 1:nfiles
%      if avgG(i) > 1.05*meanavgG | avgG(i) < 0.95*meanavgG
%         avgG(i) = meanavgG;
%      end
%  end
% 
%  for i = 1:nfiles
%      if avgB(i) > 1.05*meanavgB | avgB(i) < 0.95*meanavgB
%          avgB(i) = meanavgB;
%      end
%  end

%% FILTERING ENDS

xaxis = 1:nfiles;

%% DISPLAY RAW MEAN DATA FOR SILICON

figure(1)
hold on
plot(xaxis,meanR_0nm,'r')
plot(xaxis,meanG_0nm,'g')
plot(xaxis,meanB_0nm,'b')
xlabel('# images')
ylabel('not normalized intensity')
title('raw mean data for silicon')
subfolder_name=char(append('/CHIP',Chip_no,'_',Exposure,'ms_',FPS,'FPS',ImgCount,'_IMG',Camera,'/CHIP_',Chip_no,'_',Exposure,'ms_',FPS,'FPS',ImgCount,'_IMG',Camera,'fig1'));
saveas(figure(1),fullfile(projectdir, subfolder_name));

%% FIX SAVING AT THE FOLDER TOMORROW

figure(2)
hold on
plot(xaxis,meanR_Xnm,'r')
plot(xaxis,meanG_Xnm,'g')
plot(xaxis,meanB_Xnm,'b')
xlabel('# images')
ylabel('not normalized intensity')
title('raw mean data for silicon oxide')
subfolder_name=char(append('/CHIP',Chip_no,'_',Exposure,'ms_',FPS,'FPS',ImgCount,'_IMG',Camera,'/CHIP_',Chip_no,'_',Exposure,'ms_',FPS,'FPS',ImgCount,'_IMG',Camera,'fig2'));
saveas(figure(2),fullfile(projectdir, subfolder_name));

%% Ask user to input range for cutoff values

pause(1)
HCutoffR = input("Enter red high cutoff value: ");
LCutoffR = input("Enter red low cutoff value: ");
HCutoffG = input("Enter green high cutoff value: ");
LCutoffG = input("Enter green low cutoff value: ");
HCutoffB = input("Enter blue high cutoff value: ");
LCutoffB = input("Enter blue low cutoff value: ");
%%

R_bigger = double(meanR_0nm>HCutoffR);
R_smaller = double(meanR_0nm<LCutoffR);
G_bigger = double(meanG_0nm>HCutoffG);
G_smaller = double(meanG_0nm<LCutoffG);
B_bigger = double(meanB_0nm>HCutoffB);
B_smaller = double(meanB_0nm<LCutoffB);
%%
BigSmall = [R_bigger;R_smaller;G_bigger;G_smaller;B_bigger;B_smaller];
%%
Elims = double(find(mean(BigSmall) ~= 0));
%% DELETE VARIABLES THAT SHOULD BE ELIMINATED

filtmeanR_0nm = meanR_0nm;
filtmeanG_0nm = meanG_0nm;
filtmeanB_0nm = meanB_0nm;

filtmeanR_Xnm = meanR_Xnm;
filtmeanG_Xnm = meanG_Xnm;
filtmeanB_Xnm = meanB_Xnm;

filtmeanR_0nm(Elims) = [];
filtmeanG_0nm(Elims) = [];
filtmeanB_0nm(Elims) = [];

filtmeanR_Xnm(Elims) = [];
filtmeanG_Xnm(Elims) = [];
filtmeanB_Xnm(Elims) = [];

nonzeros_R_0nm(:,Elims) = [];
nonzeros_G_0nm(:,Elims) = [];
nonzeros_B_0nm(:,Elims) = [];

nonzeros_R_Xnm(:,Elims) = [];
nonzeros_G_Xnm(:,Elims) = [];
nonzeros_B_Xnm(:,Elims) = [];
%%

thin_ring_R_clean = thin_ring_R;
thin_ring_G_clean = thin_ring_G;
thin_ring_B_clean = thin_ring_B;

thin_ring_R_clean(:,:,Elims) = [];
thin_ring_G_clean(:,:,Elims) = [];
thin_ring_B_clean(:,:,Elims) = [];

%% RAW INTENSITY WITH ELIMINATION LINES

figure(3)
hold on
plot(xaxis,meanR_0nm,'r')
plot(xaxis,meanG_0nm,'g')
plot(xaxis,meanB_0nm,'b')
xline(Elims)
% title('FPS =')
xlabel('# images')
ylabel('not normalized intensity')
title('Silicon, vertical lines are to be eliminated')
subfolder_name=char(append('/CHIP',Chip_no,'_',Exposure,'ms_',FPS,'FPS',ImgCount,'_IMG',Camera,'/CHIP_',Chip_no,'_',Exposure,'ms_',FPS,'FPS',ImgCount,'_IMG',Camera,'fig3'));
saveas(figure(3),fullfile(projectdir, subfolder_name));
%%
figure(4)
hold on
plot(xaxis,meanR_Xnm,'r')
plot(xaxis,meanG_Xnm,'g')
plot(xaxis,meanB_Xnm,'b')
xline(Elims)
% title('FPS =')
xlabel('# images')
ylabel('not normalized intensity')
title('Silicon Oxide, vertical lines are to be eliminated')
subfolder_name=char(append('/CHIP',Chip_no,'_',Exposure,'ms_',FPS,'FPS',ImgCount,'_IMG',Camera,'/CHIP_',Chip_no,'_',Exposure,'ms_',FPS,'FPS',ImgCount,'_IMG',Camera,'fig4'));
saveas(figure(4),fullfile(projectdir, subfolder_name));
%% NOT NORMALIZED AND FILTERED INTENSITY

filtxaxis = 1:length(filtmeanB_Xnm)

figure(5)
hold on
plot(filtxaxis,filtmeanR_0nm,'r')
plot(filtxaxis,filtmeanG_0nm,'g')
plot(filtxaxis,filtmeanB_0nm,'b')
xlabel('# images')
ylabel('not normalized intensity')
title('Silicon, outliers filtered out')
subfolder_name=char(append('/CHIP',Chip_no,'_',Exposure,'ms_',FPS,'FPS',ImgCount,'_IMG',Camera,'/CHIP_',Chip_no,'_',Exposure,'ms_',FPS,'FPS',ImgCount,'_IMG',Camera,'fig5'));
saveas(figure(5),fullfile(projectdir, subfolder_name));
%%
figure(6)
hold on
plot(filtxaxis,filtmeanR_Xnm,'r')
plot(filtxaxis,filtmeanG_Xnm,'g')
plot(filtxaxis,filtmeanB_Xnm,'b')
xlabel('# images')
ylabel('not normalized intensity')
title('Silicon Oxide outliers filtered out')
subfolder_name=char(append('/CHIP',Chip_no,'_',Exposure,'ms_',FPS,'FPS',ImgCount,'_IMG',Camera,'/CHIP_',Chip_no,'_',Exposure,'ms_',FPS,'FPS',ImgCount,'_IMG',Camera,'fig6'));
saveas(figure(6),fullfile(projectdir, subfolder_name));

%% Normalize Intensity v Time and display (Normalized Filtered Average)

normfiltR_0nm = filtmeanR_0nm/max(max(filtmeanR_0nm));
normfiltG_0nm = filtmeanG_0nm/max(max(filtmeanG_0nm));
normfiltB_0nm = filtmeanB_0nm/max(max(filtmeanB_0nm));

normfiltR_Xnm = filtmeanR_Xnm/max(max(filtmeanR_Xnm));
normfiltG_Xnm = filtmeanG_Xnm/max(max(filtmeanG_Xnm));
normfiltB_Xnm = filtmeanB_Xnm/max(max(filtmeanB_Xnm));

%% NORMALIZED AND FILTERED INTENSITY

figure(7)
hold on
plot(filtxaxis,normfiltR_0nm*100,'r')
plot(filtxaxis,normfiltG_0nm*100,'g')
plot(filtxaxis,normfiltB_0nm*100,'b')
ylim([90 100])
xlim([1 length(normfiltB_0nm)])
% title('FPS =')
xlabel('# images')
ylabel('normalized intensity')
title('Silicon, outliers filtered out')
subfolder_name=char(append('/CHIP',Chip_no,'_',Exposure,'ms_',FPS,'FPS',ImgCount,'_IMG',Camera,'/CHIP_',Chip_no,'_',Exposure,'ms_',FPS,'FPS',ImgCount,'_IMG',Camera,'fig7'));
saveas(figure(7),fullfile(projectdir, subfolder_name));

%%
figure(8)
hold on
plot(filtxaxis,normfiltR_Xnm*100,'r')
plot(filtxaxis,normfiltG_Xnm*100,'g')
plot(filtxaxis,normfiltB_Xnm*100,'b')
ylim([90 100])
xlim([1 length(normfiltB_Xnm)])
% title('FPS =')
xlabel('# images')
ylabel('normalized intensity')
title('Silicon Oxide, outliers filtered out')
subfolder_name=char(append('/CHIP',Chip_no,'_',Exposure,'ms_',FPS,'FPS',ImgCount,'_IMG',Camera,'/CHIP_',Chip_no,'_',Exposure,'ms_',FPS,'FPS',ImgCount,'_IMG',Camera,'fig8'));
saveas(figure(8),fullfile(projectdir, subfolder_name));

%%
percentflucR_0nm = max(normfiltR_0nm) - min(normfiltR_0nm);
percentflucG_0nm = max(normfiltG_0nm) - min(normfiltG_0nm);
percentflucB_0nm = max(normfiltB_0nm) - min(normfiltB_0nm);

percentflucR_Xnm = max(normfiltR_Xnm) - min(normfiltR_Xnm);
percentflucG_Xnm = max(normfiltG_Xnm) - min(normfiltG_Xnm);
percentflucB_Xnm = max(normfiltB_Xnm) - min(normfiltB_Xnm);

fprintf("Percent fluctuation at Silicon at R is %f \n",percentflucR_0nm)
fprintf("Percent fluctuation at Silicon at G is %f \n",percentflucG_0nm)
fprintf("Percent fluctuation at Silicon at B is %f \n",percentflucB_0nm)

fprintf("Percent fluctuation at Silicon Oxide at R is %f \n",percentflucR_Xnm)
fprintf("Percent fluctuation at Silicon Oxide at G is %f \n",percentflucG_Xnm)
fprintf("Percent fluctuation at Silicon Oxide at B is %f \n",percentflucB_Xnm)
%%

num = length(normfiltB_0nm);
stepsize = input("Enter step size: ");

%% NOT NORMALIZED AND FILTERED TEMPORAL AVERAGE WITH STEP SIZE
for i = 1:num/stepsize
    tempavgR_0nm(i) = mean(filtmeanR_0nm((stepsize*(i-1))+1:stepsize*i));
    tempavgG_0nm(i) = mean(filtmeanG_0nm((stepsize*(i-1))+1:stepsize*i));
    tempavgB_0nm(i) = mean(filtmeanB_0nm((stepsize*(i-1))+1:stepsize*i));
end
%%
for i = 1:num/stepsize
    tempavgR_Xnm(i) = mean(filtmeanR_Xnm((stepsize*(i-1))+1:stepsize*i));
    tempavgG_Xnm(i) = mean(filtmeanG_Xnm((stepsize*(i-1))+1:stepsize*i));
    tempavgB_Xnm(i) = mean(filtmeanB_Xnm((stepsize*(i-1))+1:stepsize*i));
end

%% FOR THICKNESS CALCULATION

for i = 1:num/stepsize
    tempavgR_0nm_th(:,:,i) = uint16(sum(thin_ring_R_clean(:,:,(stepsize*(i-1))+1:stepsize*i),3) ./ stepsize);
    tempavgG_0nm_th(:,:,i) = uint16(sum(thin_ring_G_clean(:,:,(stepsize*(i-1))+1:stepsize*i),3) ./ stepsize);
    tempavgB_0nm_th(:,:,i) = uint16(sum(thin_ring_B_clean(:,:,(stepsize*(i-1))+1:stepsize*i),3) ./ stepsize);
end
%%
for i = 1:num/stepsize
    tempavgR_Xnm_th(:,:,i) = uint16(sum(cropped_im_R_Xnm(:,:,(stepsize*(i-1))+1:stepsize*i),3) ./ stepsize);
    tempavgG_Xnm_th(:,:,i) = uint16(sum(cropped_im_G_Xnm(:,:,(stepsize*(i-1))+1:stepsize*i),3) ./ stepsize);
    tempavgB_Xnm_th(:,:,i) = uint16(sum(cropped_im_B_Xnm(:,:,(stepsize*(i-1))+1:stepsize*i),3) ./ stepsize);
end

% nonzeros_0nm

%%
tempavgxaxis = 1:length(tempavgB_0nm);

figure(9)
hold on
plot(tempavgxaxis,tempavgR_0nm,'r')
plot(tempavgxaxis,tempavgG_0nm,'g')
plot(tempavgxaxis,tempavgB_0nm,'b')
ylim([0 max(tempavgB_0nm)])
xlim([1 length(tempavgB_0nm)])
xlabel('# images')
ylabel('not normalized intensity')
title('Silicon, outliers filtered out')
subfolder_name=char(append('/CHIP',Chip_no,'_',Exposure,'ms_',FPS,'FPS',ImgCount,'_IMG',Camera,'/CHIP_',Chip_no,'_',Exposure,'ms_',FPS,'FPS',ImgCount,'_IMG',Camera,'fig9'));
saveas(figure(9),fullfile(projectdir, subfolder_name));

%%

figure(10)
hold on
plot(tempavgxaxis,tempavgR_Xnm,'r')
plot(tempavgxaxis,tempavgG_Xnm,'g')
plot(tempavgxaxis,tempavgB_Xnm,'b')
ylim([0 max(tempavgB_Xnm)])
xlim([1 length(tempavgB_Xnm)])
xlabel('# images')
ylabel('not normalized intensity')
title('Silicon Oxide, outliers filtered out')
subfolder_name=char(append('/CHIP',Chip_no,'_',Exposure,'ms_',FPS,'FPS',ImgCount,'_IMG',Camera,'/CHIP_',Chip_no,'_',Exposure,'ms_',FPS,'FPS',ImgCount,'_IMG',Camera,'fig10'));
saveas(figure(10),fullfile(projectdir, subfolder_name));

%% NORMALIZED AND FILTERED TEMPORAL AVERAGE WITH STEPSIZE

for i = 1:num/stepsize
    tempavgRnorm_0nm(i) = 100*mean(normfiltR_0nm((stepsize*(i-1))+1:stepsize*i));
    tempavgGnorm_0nm(i) = 100*mean(normfiltG_0nm((stepsize*(i-1))+1:stepsize*i));
    tempavgBnorm_0nm(i) = 100*mean(normfiltB_0nm((stepsize*(i-1))+1:stepsize*i));
end

for i = 1:num/stepsize
    tempavgRnorm_Xnm(i) = 100*mean(normfiltR_Xnm((stepsize*(i-1))+1:stepsize*i));
    tempavgGnorm_Xnm(i) = 100*mean(normfiltG_Xnm((stepsize*(i-1))+1:stepsize*i));
    tempavgBnorm_Xnm(i) = 100*mean(normfiltB_Xnm((stepsize*(i-1))+1:stepsize*i));
end

%%
figure(11)
hold on
plot(tempavgxaxis,tempavgRnorm_0nm,'r')
plot(tempavgxaxis,tempavgGnorm_0nm,'g')
plot(tempavgxaxis,tempavgBnorm_0nm,'b')
ylim([90 100])
xlim([1 length(tempavgB_0nm)])
xlabel('# images')
ylabel('normalized intensity')
title('Silicon, outliers filtered out')
subfolder_name=char(append('/CHIP',Chip_no,'_',Exposure,'ms_',FPS,'FPS',ImgCount,'_IMG',Camera,'/CHIP_',Chip_no,'_',Exposure,'ms_',FPS,'FPS',ImgCount,'_IMG',Camera,'fig11'));
saveas(figure(11),fullfile(projectdir, subfolder_name));
%%
figure(12)
hold on
plot(tempavgxaxis,tempavgRnorm_Xnm,'r')
plot(tempavgxaxis,tempavgGnorm_Xnm,'g')
plot(tempavgxaxis,tempavgBnorm_Xnm,'b')
ylim([90 100])
xlim([1 length(tempavgB_Xnm)])
xlabel('# images')
ylabel('normalized intensity')
title('Silicon Oxide, outliers filtered out')
subfolder_name=char(append('/CHIP',Chip_no,'_',Exposure,'ms_',FPS,'FPS',ImgCount,'_IMG',Camera,'/CHIP_',Chip_no,'_',Exposure,'ms_',FPS,'FPS',ImgCount,'_IMG',Camera,'fig12'));
saveas(figure(12),fullfile(projectdir, subfolder_name));
%%


for i = 1:size(tempavgB_Xnm_th,3)
nonzeros_R_0nm_th(:,i) = nonzeros(tempavgR_0nm_th(:,:,i));
nonzeros_G_0nm_th(:,i) = nonzeros(tempavgG_0nm_th(:,:,i));
nonzeros_B_0nm_th(:,i) = nonzeros(tempavgB_0nm_th(:,:,i));

nonzeros_R_Xnm_th(:,i) = nonzeros(tempavgR_Xnm_th(:,:,i));
nonzeros_G_Xnm_th(:,i) = nonzeros(tempavgG_Xnm_th(:,:,i));
nonzeros_B_Xnm_th(:,i) = nonzeros(tempavgB_Xnm_th(:,:,i));
end


%%

if numel(nonzeros_R_0nm_th) > numel(nonzeros_R_Xnm_th)
    nonzeros_R_0nm_th = nonzeros_R_0nm_th(1:rowsnonzerosXnm,:);
elseif numel(nonzeros_R_0nm_th) < numel(nonzeros_R_Xnm_th)
    nonzeros_R_Xnm_th = nonzeros_R_Xnm_th(1:rowsnonzeros0nm,:); 
end
%%
if numel(nonzeros_G_0nm_th) > numel(nonzeros_G_Xnm_th)
    nonzeros_G_0nm_th = nonzeros_G_0nm_th(1:rowsnonzerosXnm,:);
elseif numel(nonzeros_G_0nm_th) < numel(nonzeros_G_Xnm_th)
    nonzeros_G_Xnm_th = nonzeros_G_Xnm_th(1:rowsnonzeros0nm,:); 
end

if numel(nonzeros_B_0nm_th) > numel(nonzeros_B_Xnm_th)
    nonzeros_B_0nm_th = nonzeros_B_0nm_th(1:rowsnonzerosXnm,:);
elseif numel(nonzeros_B_0nm_th) < numel(nonzeros_B_Xnm_th)
    nonzeros_B_Xnm_th = nonzeros_B_Xnm_th(1:rowsnonzeros0nm,:); 
end




%% Intensity_in calculation (calibration) I_in = I_out_Si / Reflectance

I_inR = round(double(nonzeros_R_0nm_th) ./ Ref_at_0nm_R);
I_inG = round(double(nonzeros_G_0nm_th) ./ Ref_at_0nm_G);
I_inB = round(double(nonzeros_B_0nm_th) ./ Ref_at_0nm_B);

%% Reflectance from oxide calculation: Reflectance_Oxide = I_Out_SiO2 / I_in

RefRed_at_XnmO2 = (double(nonzeros_R_Xnm_th)./I_inR);
RefGre_at_XnmO2 = (double(nonzeros_G_Xnm_th)./I_inG);
RefBlu_at_XnmO2 = (double(nonzeros_B_Xnm_th)./I_inB);

%% Take the average of Reflectance values in R,G,B only pixels

MeanRefRed_at_XnmO2 = mean(mean(RefRed_at_XnmO2));
MeanRefGre_at_XnmO2 = mean(mean(RefGre_at_XnmO2));
MeanRefBlu_at_XnmO2 = mean(mean(RefBlu_at_XnmO2));

%% Use reftocurve function to estimate thickness from these reflectance values

esti_L = reftocurve_lsqr(MeanRefRed_at_XnmO2,MeanRefGre_at_XnmO2,MeanRefBlu_at_XnmO2);

%% Display estimated thicknesses and fitted thickness curves for RGB

figure(13)
hold on

%% Plot calculated reflectance values at RGB

plot(cw_r,MeanRefRed_at_XnmO2,'r.','MarkerSize',30) 
plot(cw_g,MeanRefGre_at_XnmO2,'g.','MarkerSize',30)
plot(cw_b,MeanRefBlu_at_XnmO2,'b.','MarkerSize',30) 

%% Plot fitted reflectance curves at estimated thickness for RGB

plot(lambda,Gamma_B(:,valtoindex_L(abs(esti_L))),'b--','LineWidth',2)
plot(lambda,Gamma_G(:,valtoindex_L(abs(esti_L))),'g--','LineWidth',2)
plot(lambda,Gamma_R(:,valtoindex_L(abs(esti_L))),'r--','LineWidth',2)
% plot(lambda,Gamma(:,valtoindex_L(abs(esti_L))),'c--','LineWidth',2)

%% Plot simulated reflectance curve at the theoretical thickness for RGB

plot(lambda,Gamma_B(:,valtoindex_L(theoric_L)),'b--','LineWidth',2)
plot(lambda,Gamma_G(:,valtoindex_L(theoric_L)),'g--','LineWidth',2)
plot(lambda,Gamma_R(:,valtoindex_L(theoric_L)),'r--','LineWidth',2)
% plot(lambda,Gamma(:,valtoindex_L(theoric_L)),'c--','LineWidth',2)

xlabel('lambda (\mum)')
ylabel('Reflectance')
xlim([0.4 0.68])
ylim([0 1])
% legend1 = sprintf('B Curve at L = %0.2f nm', 1000*esti_L);
% legend2 = sprintf('G Curve at L = %0.2f nm', 1000*esti_L);
% legend3 = sprintf('R Curve at L = %0.2f nm', 1000*esti_L);
% legend4 = sprintf('B Curve at L = %0.2f nm', 1000*theoric_L);
% legend5 = sprintf('G Curve at L = %0.2f nm', 1000*theoric_L);
% legend6 = sprintf('R Curve at L = %0.2f nm', 1000*theoric_L);

legend1 = sprintf('B Curves', 1000*esti_L);
legend2 = sprintf('G Curves', 1000*esti_L);
legend3 = sprintf('R Curves', 1000*esti_L);

legend('Ref at R','Ref at G','Ref at B','B Curves', 'G Curves', 'R Curves','location','bestoutside')
title("Actual Thickness is " + theoric_L*1000 + " nm, Estimated Thickness is " + esti_L*1000 + " nm")

subfolder_name=char(append('/CHIP',Chip_no,'_',Exposure,'ms_',FPS,'FPS',ImgCount,'_IMG',Camera,'/CHIP_',Chip_no,'_',Exposure,'ms_',FPS,'FPS',ImgCount,'_IMG',Camera,'fig13'));
saveas(figure(13),fullfile(projectdir, subfolder_name));

if isnumeric(theoric_L)
percent_error = 100*((abs(esti_L-theoric_L))/theoric_L);
fprintf('Actual thickness is %f nm \n',theoric_L*1000)
fprintf('Estimated thickness is %f nm \n',esti_L*1000)
fprintf("Percent error is %f \", percent_error)
% elseif theoric_L == 'UNKNOWN'
% fprintf('Estimated thickness is %f nm \n',esti_L*1000)
% disp('Percent Error is not known')
end

%% SUPPLEMENTARY IMAGES AND FIGURES

%% CODE TO SHOW MANUALLY SELECTED OXIDE OF FIRST IMAGE OF SELECTION (UN/COMMENT STARTING HERE)

figure(14)
imagesc(cropped_im_R_Xnm(:,:,1))
colormap
colorbar
caxis([min(min(nonzeros(cropped_im_R_Xnm))) max(max(nonzeros(cropped_im_R_Xnm)))]);
title('Cropped R Xnm ROI')
subfolder_name=char(append('/CHIP',Chip_no,'_',Exposure,'ms_',FPS,'FPS',ImgCount,'_IMG',Camera,'/CHIP_',Chip_no,'_',Exposure,'ms_',FPS,'FPS',ImgCount,'_IMG',Camera,'fig14'));
saveas(figure(14),fullfile(projectdir, subfolder_name));

figure(15)
imagesc(cropped_im_G_Xnm(:,:,1))
colormap
colorbar
caxis([min(min(nonzeros(cropped_im_G_Xnm))) max(max(nonzeros(cropped_im_G_Xnm)))]);
title('Cropped G Xnm ROI')
subfolder_name=char(append('/CHIP',Chip_no,'_',Exposure,'ms_',FPS,'FPS',ImgCount,'_IMG',Camera,'/CHIP_',Chip_no,'_',Exposure,'ms_',FPS,'FPS',ImgCount,'_IMG',Camera,'fig15'));
saveas(figure(15),fullfile(projectdir, subfolder_name));

figure(16)
imagesc(cropped_im_B_Xnm(:,:,1))
colormap
colorbar
caxis([min(min(nonzeros(cropped_im_B_Xnm))) max(max(nonzeros(cropped_im_B_Xnm)))]);
title('Cropped B Xnm ROI')
subfolder_name=char(append('/CHIP',Chip_no,'_',Exposure,'ms_',FPS,'FPS',ImgCount,'_IMG',Camera,'/CHIP_',Chip_no,'_',Exposure,'ms_',FPS,'FPS',ImgCount,'_IMG',Camera,'fig16'));
saveas(figure(16),fullfile(projectdir, subfolder_name));

%% UN/COMMENT ENDS

%% CODE TO SHOW RGB RINGS OF FIRST IMAGE OF SELECTION

figure(17)
imagesc(thin_ring_R(:,:,1))
colormap
colorbar
caxis([min(min(nonzeros(thin_ring_R))) max(max(nonzeros(thin_ring_R)))]);
title('R Ring ROI')
subfolder_name=char(append('/CHIP',Chip_no,'_',Exposure,'ms_',FPS,'FPS',ImgCount,'_IMG',Camera,'/CHIP_',Chip_no,'_',Exposure,'ms_',FPS,'FPS',ImgCount,'_IMG',Camera,'fig17'));
saveas(figure(17),fullfile(projectdir, subfolder_name));

figure(18)
imagesc(thin_ring_G(:,:,1))
colormap
colorbar
caxis([min(min(nonzeros(thin_ring_G))) max(max(nonzeros(thin_ring_G)))]);
title('G Ring ROI')
subfolder_name=char(append('/CHIP',Chip_no,'_',Exposure,'ms_',FPS,'FPS',ImgCount,'_IMG',Camera,'/CHIP_',Chip_no,'_',Exposure,'ms_',FPS,'FPS',ImgCount,'_IMG',Camera,'fig18'));
saveas(figure(18),fullfile(projectdir, subfolder_name));

figure(19)
imagesc(thin_ring_B(:,:,1))
colormap
colorbar
caxis([double(min(min(nonzeros(thin_ring_B)))) double(max(max(nonzeros(thin_ring_B))))]);
title('B Ring ROI')
subfolder_name=char(append('/CHIP',Chip_no,'_',Exposure,'ms_',FPS,'FPS',ImgCount,'_IMG',Camera,'/CHIP_',Chip_no,'_',Exposure,'ms_',FPS,'FPS',ImgCount,'_IMG',Camera,'fig19'));
saveas(figure(19),fullfile(projectdir, subfolder_name));

%% UN/COMMENT ENDS



%%
toc