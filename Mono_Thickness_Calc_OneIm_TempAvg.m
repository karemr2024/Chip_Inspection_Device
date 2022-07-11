%% START

tic
clc; clearvars; close all;
fprintf('Beginning to run %s.m ...\n', mfilename);

%% Below are prompts for user to input unknown SiO2 nm images in .mat format

[R_Tiff_Name,R_Tiff_Path] = uigetfile('*.mat','Red Stack'); %Import Unknown Thickness: Red Image
if isequal(R_Tiff_Name,0)
   disp('User selected Cancel');
else
   disp(['User selected ', fullfile(R_Tiff_Name,R_Tiff_Path)]);
   R_Stack_pre = load(strcat(R_Tiff_Path, R_Tiff_Name));
end

[G_Tiff_Name,G_Tiff_Path] = uigetfile('*.mat','Green Stack'); %Import Unknown Thickness: Green Image
if isequal(G_Tiff_Name,0)
   disp('User selected Cancel');
else
   disp(['User selected ', fullfile(G_Tiff_Name,G_Tiff_Path)]);
   G_Stack_pre = load(strcat(G_Tiff_Path, G_Tiff_Name));
end

[B_Tiff_Name,B_Tiff_Path] = uigetfile('*.mat','Blue Stack'); %Import Unknown Thickness: Blue Image
if isequal(B_Tiff_Name,0)
   disp('User selected Cancel');
else
   disp(['User selected ', fullfile(B_Tiff_Name,B_Tiff_Path)]);
   B_Stack_pre = load(strcat(B_Tiff_Path, B_Tiff_Name));
end

%% Prompt user about the thickness of the chip

load("Simulation_Data.mat")
prompt_NorS = input('Do you know the thickness of your chip? (Y/N)\n', "s");
if prompt_NorS == 'Y'
theoric_L = (input('Enter the thickness of your chip in nm: \n'))/1000;
else
theoric_L = 'UNKNOWN';
end

%% Convert RGB stacks from struct to cell

R_Stack = struct2cell(R_Stack_pre);
G_Stack = struct2cell(G_Stack_pre);
B_Stack = struct2cell(B_Stack_pre);

%% Separate images and timestamps

Rims = R_Stack(1:2:end);
Rtimes =R_Stack(2:2:end);
Gims = G_Stack(1:2:end);
Gtimes = G_Stack(2:2:end);
Bims = B_Stack(1:2:end);
Btimes = B_Stack(2:2:end);

%% Get size

[rows, cols] = size(cell2mat(Rims(1)));
nfiles = length(Rims);

%% Convert RGB stacks to matrices

tiff_stack_R = uint16(zeros(rows,cols,nfiles));
tiff_stack_G = uint16(zeros(rows,cols,nfiles));
tiff_stack_B = uint16(zeros(rows,cols,nfiles));
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

%% Display ROI selected spot image and select ROI for outer Silicon ring

[cropped_im_R_0nm_out(:,:,1),x_crop_0nm_out,y_crop_0nm_out] = AreaSelection_Circle_SiOut(R_Crop_Bothnm(:,:,1));
cropped_im_G_0nm_out(:,:,1) = AreaSelection_Circle_Mod(G_Crop_Bothnm(:,:,1),x_crop_0nm_out,y_crop_0nm_out);
cropped_im_B_0nm_out(:,:,1) = AreaSelection_Circle_Mod(B_Crop_Bothnm(:,:,1),x_crop_0nm_out,y_crop_0nm_out);

for i = 2:nfiles
    cropped_im_R_0nm_out(:,:,i) = AreaSelection_Circle_Mod(R_Crop_Bothnm(:,:,i),x_crop_0nm_out,y_crop_0nm_out);
    cropped_im_G_0nm_out(:,:,i) = AreaSelection_Circle_Mod(G_Crop_Bothnm(:,:,i),x_crop_0nm_out,y_crop_0nm_out);
    cropped_im_B_0nm_out(:,:,i) = AreaSelection_Circle_Mod(B_Crop_Bothnm(:,:,i),x_crop_0nm_out,y_crop_0nm_out);
end

%% Display ROI selected outer ring image and select ROI for inner Silicon ring

[cropped_im_R_0nm_in(:,:,1),x_crop_0nm_in,y_crop_0nm_in] = AreaSelection_Circle_SiIn(cropped_im_R_0nm_out(:,:,1),x_crop_0nm_out,y_crop_0nm_out);
cropped_im_G_0nm_in(:,:,1) = AreaSelection_Circle_Mod_SiIn(cropped_im_G_0nm_out(:,:,1),x_crop_0nm_in,y_crop_0nm_in);
cropped_im_B_0nm_in(:,:,1) = AreaSelection_Circle_Mod_SiIn(cropped_im_B_0nm_out(:,:,1),x_crop_0nm_in,y_crop_0nm_in);

for i = 2:nfiles
    cropped_im_R_0nm_in(:,:,i) = AreaSelection_Circle_Mod_SiIn(cropped_im_R_0nm_out(:,:,i),x_crop_0nm_in,y_crop_0nm_in);
    cropped_im_G_0nm_in(:,:,i) = AreaSelection_Circle_Mod_SiIn(cropped_im_G_0nm_out(:,:,i),x_crop_0nm_in,y_crop_0nm_in);
    cropped_im_B_0nm_in(:,:,i) = AreaSelection_Circle_Mod_SiIn(cropped_im_B_0nm_out(:,:,i),x_crop_0nm_in,y_crop_0nm_in);
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

%% INTENSITY v TIME GRAPHS AT SILICON ONLY (THIN RING)

%% Get the mean of the silicon oxide values

meanavgR = mean(nonzeros_R_0nm);
meanavgG = mean(nonzeros_G_0nm);
meanavgB = mean(nonzeros_B_0nm);

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

figure
hold on
plot(xaxis,meanavgR,'r')
plot(xaxis,meanavgG,'g')
plot(xaxis,meanavgB,'b')
% title('FPS =')
xlabel('# images')
ylabel('not normalized intensity')

%% Ask user to input range for temporal average
pause(1)
firstim = input('Begin Temporal Average = ');
pause(1)
lastim =  input('End Temporal Average = ');

%% Choose sample for inputted range FIX HERE!!!!!!!!!!!!!
close(figure)

% sample_R = zeros(rows,cols,lastim-firstim+1);
% sample_G = zeros(rows,cols,lastim-firstim+1);
% sample_B = zeros(rows,cols,lastim-firstim+1); %%OLD VERSION

% for i = 1:lastim-firstim+1
% sample_R(:,:,i) = cell2mat(Rims(firstim-1+i));
% sample_G(:,:,i) = cell2mat(Gims(firstim-1+i)); 
% sample_B(:,:,i) = cell2mat(Bims(firstim-1+i)); 
% end 

for i = 1:lastim-firstim+1
sample_R(:,:,i) = tiff_stack_R(firstim-1+i); %NEW VERSION UPDATE TOMORROW
sample_G(:,:,i) = tiff_stack_G(firstim-1+i); 
sample_B(:,:,i) = tiff_stack_B(firstim-1+i); 
end 

%% Normalize Intensity v Time and display

normmeanavgsampleR = meanavgR(firstim:lastim);
normmeanavgsampleG = meanavgG(firstim:lastim);
normmeanavgsampleB = meanavgB(firstim:lastim);

normmeanavgR = normmeanavgsampleR/max(max(normmeanavgsampleR));
normmeanavgG = normmeanavgsampleG/max(max(normmeanavgsampleG));
normmeanavgB = normmeanavgsampleB/max(max(normmeanavgsampleB));
%%
normxaxis = firstim:lastim;

figure
hold on
plot(normxaxis,normmeanavgR*100,'r')
plot(normxaxis,normmeanavgG*100,'g')
plot(normxaxis,normmeanavgB*100,'b')
ylim([90 100])
xlim([firstim lastim])
% title('FPS =')
xlabel('# images')
ylabel('normalized intensity')

%%

num = length(normmeanavgB);
stepsize = input("Enter step size: ")

for i = 1:num/stepsize
    tempavgR(i) = 100*mean(normmeanavgR(i*10-9:stepsize*i))
    tempavgG(i) = 100*mean(normmeanavgG(i*10-9:stepsize*i))
    tempavgB(i) = 100*mean(normmeanavgB(i*10-9:stepsize*i))
end

tempavgxaxis = firstim:stepsize:lastim;
%%
figure
hold on
plot(tempavgxaxis,tempavgR,'r')
plot(tempavgxaxis,tempavgG,'g')
plot(tempavgxaxis,tempavgB,'b')
ylim([90 100])
xlim([firstim lastim])
xlabel('# images')
ylabel('normalized intensity')

%%
percentflucR = max(normmeanavgR) - min(normmeanavgR)
percentflucG = max(normmeanavgG) - min(normmeanavgG)
percentflucB = max(normmeanavgB) - min(normmeanavgB)

fprintf("Percent fluctuation at R is %f \n",percentflucR)
fprintf("Percent fluctuation at G is %f \n",percentflucG)
fprintf("Percent fluctuation at B is %f \n",percentflucB)
%% Make the number of elements in silicon and oxide equal

if numel(nonzeros_R_0nm) > numel(nonzeros_R_Xnm)
    nonzeros_R_0nm = nonzeros_R_0nm(1:rowsnonzerosXnm,:);
elseif numel(nonzeros_R_0nm) < numel(nonzeros_R_Xnm)
    nonzeros_R_Xnm = nonzeros_R_Xnm(1:rowsnonzeros0nm,:); 
end

if numel(nonzeros_G_0nm) > numel(nonzeros_G_Xnm)
    nonzeros_G_0nm = nonzeros_G_0nm(1:rowsnonzerosXnm,:);
elseif numel(nonzeros_G_0nm) < numel(nonzeros_G_Xnm)
    nonzeros_G_Xnm = nonzeros_G_Xnm(1:rowsnonzeros0nm,:); 
end

if numel(nonzeros_B_0nm) > numel(nonzeros_B_Xnm)
    nonzeros_B_0nm = nonzeros_B_0nm(1:rowsnonzerosXnm,:);
elseif numel(nonzeros_B_0nm) < numel(nonzeros_B_Xnm)
    nonzeros_B_Xnm = nonzeros_B_Xnm(1:rowsnonzeros0nm,:); 
end

%% Intensity_in calculation (calibration) I_in = I_out_Si / Reflectance

I_inR = round(double(nonzeros_R_0nm) ./ Ref_at_0nm_R);
I_inG = round(double(nonzeros_G_0nm) ./ Ref_at_0nm_G);
I_inB = round(double(nonzeros_B_0nm) ./ Ref_at_0nm_B);

%% Reflectance from oxide calculation: Reflectance_Oxide = I_Out_SiO2 / I_in

RefRed_at_XnmO2 = (double(nonzeros_R_Xnm)./I_inR);
RefGre_at_XnmO2 = (double(nonzeros_G_Xnm)./I_inG);
RefBlu_at_XnmO2 = (double(nonzeros_B_Xnm)./I_inB);

%% Take the average of Reflectance values in R,G,B only pixels

MeanRefRed_at_XnmO2 = mean(mean(RefRed_at_XnmO2));
MeanRefGre_at_XnmO2 = mean(mean(RefGre_at_XnmO2));
MeanRefBlu_at_XnmO2 = mean(mean(RefBlu_at_XnmO2));

%% Use reftocurve function to estimate thickness from these reflectance values

esti_L = reftocurve_lsqr(MeanRefRed_at_XnmO2,MeanRefGre_at_XnmO2,MeanRefBlu_at_XnmO2);

%% Display estimated thicknesses and fitted thickness curves for RGB

figure
hold on

%% Plot calculated reflectance values at RGB

plot(cw_r,MeanRefRed_at_XnmO2,'r.','MarkerSize',30) 
plot(cw_g,MeanRefGre_at_XnmO2,'g.','MarkerSize',30)
plot(cw_b,MeanRefBlu_at_XnmO2,'b.','MarkerSize',30) 

%% Plot fitted reflectance curves at estimated thickness for RGB

plot(lambda,Gamma_B(:,valtoindex_L(abs(esti_L))),'r--','LineWidth',2)
plot(lambda,Gamma_G(:,valtoindex_L(abs(esti_L))),'g--','LineWidth',2)
plot(lambda,Gamma_R(:,valtoindex_L(abs(esti_L))),'b--','LineWidth',2)
% plot(lambda,Gamma(:,valtoindex_L(abs(esti_L))),'b--','LineWidth',2)

%% Plot simulated reflectance curve at the theoretical thickness for RGB

plot(lambda,Gamma_B(:,valtoindex_L(theoric_L)),'c--','LineWidth',2)
plot(lambda,Gamma_G(:,valtoindex_L(theoric_L)),'c--','LineWidth',2)
plot(lambda,Gamma_R(:,valtoindex_L(theoric_L)),'c--','LineWidth',2)
% plot(lambda,Gamma(:,valtoindex_L(theoric_L)),'c--','LineWidth',2)

xlabel('lambda (\mum)')
ylabel('Reflectance')
xlim([0.4 0.68])
ylim([0 1])
legend1 = sprintf('L = %0.2f nm', 1000*esti_L);
legend2 = sprintf('L = %0.2f nm', 1000*theoric_L);
legend(legend1, legend2,'Ref at R','Ref at G','Ref at B','location','bestoutside')
title("Actual Thickness is " + theoric_L*1000 + " nm, Estimated Thickness is " + esti_L*1000 + " nm")

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

figure
imagesc(cropped_im_R_Xnm(:,:,1))
colormap
colorbar
caxis([min(min(nonzeros(cropped_im_R_Xnm))) max(max(nonzeros(cropped_im_R_Xnm)))]);
title('Cropped R Xnm ROI')

figure
imagesc(cropped_im_G_Xnm(:,:,1))
colormap
colorbar
caxis([min(min(nonzeros(cropped_im_G_Xnm))) max(max(nonzeros(cropped_im_G_Xnm)))]);
title('Cropped G Xnm ROI')

figure
imagesc(cropped_im_B_Xnm(:,:,1))
colormap
colorbar
caxis([min(min(nonzeros(cropped_im_B_Xnm))) max(max(nonzeros(cropped_im_B_Xnm)))]);
title('Cropped B Xnm ROI')

%% UN/COMMENT ENDS

%% CODE TO SHOW RGB RINGS OF FIRST IMAGE OF SELECTION

figure
imagesc(thin_ring_R(:,:,1))
colormap
colorbar
caxis([min(min(nonzeros(thin_ring_R))) max(max(nonzeros(thin_ring_R)))]);
title('R Ring ROI')
figure
imagesc(thin_ring_G(:,:,1))
colormap
colorbar
caxis([min(min(nonzeros(thin_ring_G))) max(max(nonzeros(thin_ring_G)))]);
title('G Ring ROI')
figure
imagesc(thin_ring_B(:,:,1))
colormap
colorbar
caxis([double(min(min(nonzeros(thin_ring_B)))) double(max(max(nonzeros(thin_ring_B))))]);
title('B Ring ROI')

%% UN/COMMENT ENDS

toc