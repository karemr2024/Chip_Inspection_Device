tic
% Get list of all BMP files in this directory
% DIR returns as a structure array.  You will need to use () and . to get
% the file names.
% Only works for files taken within the same day
clearvars; clc; close all;
%START NEW CODE

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

tic

 clearvars -except R_Stack_pre B_Stack_pre G_Stack_pre 
%% INPUT EXPERIMENT DETAILS EVERY TIME

FPS = 55
BufferMode = 'First Only'
StartEventType = 'Frame Start' 
EndEventType = 'Exposure End'
Camera = 'BFS-U3-17S7M-C' 
LED = 'RGB'

%%

R_Stack = struct2cell(R_Stack_pre);
G_Stack = struct2cell(G_Stack_pre);
B_Stack = struct2cell(B_Stack_pre);

Rims = R_Stack(1:2:end);
Rtimes =R_Stack(2:2:end);
Gims = G_Stack(1:2:end);
Gtimes = G_Stack(2:2:end);
Bims = B_Stack(1:2:end);
Btimes = B_Stack(2:2:end);
%%
% tiff_stack_R = uint16(rand(rows,cols,nfiles));
% tiff_stack_G = uint16(rand(rows,cols,nfiles));
% tiff_stack_B = uint16(rand(rows,cols,nfiles));

%%
for i = 1:nfiles
   tiff_stack_R(:,:,i) = cat(3,cell2mat(Rims(i)));
   tiff_stack_G(:,:,i) = cat(3,cell2mat(Gims(i)));
   tiff_stack_B(:,:,i) = cat(3,cell2mat(Bims(i)));
end

%%

%END NEW CODE
% imagefiles = dir('**//*.tif'); 
%
clearvars -except tiff_stack_B tiff_stack_G tiff_stack_R
[rows, cols, nfiles] = size(tiff_stack_B);
% nfiles = size(3);
%

% intensities = zeros(1,nfiles);

% Sort images according to date and time
% imagefiles = imagefiles(~[imagefiles.isdir]);
% [~,idx] = sort([imagefiles.datenum]);
% imagefiles = imagefiles(idx);
%
% [display_image_crop,x_crop_0nm,y_crop_0nm] = AreaSelection(imread(imagefiles(2).name));

[R_Crop_Bothnm(:,:,1),x_crop_Bothnm,y_crop_Bothnm] = AreaSelection_Circle(tiff_stack_R(:,:,1));
G_Crop_Bothnm(:,:,1) = AreaSelection_Circle_Mod(tiff_stack_G(:,:,1),x_crop_Bothnm,y_crop_Bothnm);
B_Crop_Bothnm(:,:,1) = AreaSelection_Circle_Mod(tiff_stack_B(:,:,1),x_crop_Bothnm,y_crop_Bothnm);

for i = 2:nfiles
    R_Crop_Bothnm(:,:,i) = AreaSelection_Circle_Mod(tiff_stack_R(:,:,i),x_crop_Bothnm,y_crop_Bothnm);
    G_Crop_Bothnm(:,:,i) = AreaSelection_Circle_Mod(tiff_stack_G(:,:,i),x_crop_Bothnm,y_crop_Bothnm);
    B_Crop_Bothnm(:,:,i) = AreaSelection_Circle_Mod(tiff_stack_B(:,:,i),x_crop_Bothnm,y_crop_Bothnm);
end
%
[cropped_im_R_0nm_out(:,:,1),x_crop_0nm_out,y_crop_0nm_out] = AreaSelection_Circle_SiOut(R_Crop_Bothnm(:,:,1));
cropped_im_G_0nm_out(:,:,1) = AreaSelection_Circle_Mod(G_Crop_Bothnm(:,:,1),x_crop_0nm_out,y_crop_0nm_out);
cropped_im_B_0nm_out(:,:,1) = AreaSelection_Circle_Mod(B_Crop_Bothnm(:,:,1),x_crop_0nm_out,y_crop_0nm_out);

for i = 2:nfiles
    cropped_im_R_0nm_out(:,:,i) = AreaSelection_Circle_Mod(R_Crop_Bothnm(:,:,i),x_crop_0nm_out,y_crop_0nm_out);
    cropped_im_G_0nm_out(:,:,i) = AreaSelection_Circle_Mod(G_Crop_Bothnm(:,:,i),x_crop_0nm_out,y_crop_0nm_out);
    cropped_im_B_0nm_out(:,:,i) = AreaSelection_Circle_Mod(B_Crop_Bothnm(:,:,i),x_crop_0nm_out,y_crop_0nm_out);
end
%
[cropped_im_R_0nm_in(:,:,1),x_crop_0nm_in,y_crop_0nm_in] = AreaSelection_Circle_SiIn(cropped_im_R_0nm_out(:,:,1),x_crop_0nm_out,y_crop_0nm_out);
cropped_im_G_0nm_in(:,:,1) = AreaSelection_Circle_Mod_SiIn(cropped_im_G_0nm_out(:,:,1),x_crop_0nm_in,y_crop_0nm_in);
cropped_im_B_0nm_in(:,:,1) = AreaSelection_Circle_Mod_SiIn(cropped_im_B_0nm_out(:,:,1),x_crop_0nm_in,y_crop_0nm_in);

for i = 2:nfiles
    cropped_im_R_0nm_in(:,:,i) = AreaSelection_Circle_Mod_SiIn(cropped_im_R_0nm_out(:,:,i),x_crop_0nm_in,y_crop_0nm_in);
    cropped_im_G_0nm_in(:,:,i) = AreaSelection_Circle_Mod_SiIn(cropped_im_G_0nm_out(:,:,i),x_crop_0nm_in,y_crop_0nm_in);
    cropped_im_B_0nm_in(:,:,i) = AreaSelection_Circle_Mod_SiIn(cropped_im_B_0nm_out(:,:,i),x_crop_0nm_in,y_crop_0nm_in);
end

%
for i = 1:nfiles
thin_ring_R(:,:,i) = uint16(cropped_im_R_0nm_out(:,:,i)) - uint16(cropped_im_R_0nm_in(:,:,i)); 
thin_ring_G(:,:,i) = uint16(cropped_im_G_0nm_out(:,:,i)) - uint16(cropped_im_G_0nm_in(:,:,i));
thin_ring_B(:,:,i) = uint16(cropped_im_B_0nm_out(:,:,i)) - uint16(cropped_im_B_0nm_in(:,:,i));
end

%

for i = 1:nfiles
    nonzeros_R(:,i) = nonzeros(thin_ring_R(:,:,i));
    nonzeros_G(:,i) = nonzeros(thin_ring_G(:,:,i));
    nonzeros_B(:,i) = nonzeros(thin_ring_B(:,:,i));
end
%
for i = 1:nfiles
    avgR(i) = mean(nonzeros_R(:,i));   
    avgG(i) = mean(nonzeros_G(:,i));   
    avgB(i) = mean(nonzeros_B(:,i)); 
end

meanavgR = mean(avgR);
meanavgG = mean(avgG);
meanavgB = mean(avgB);

%
 for i = 1:nfiles
     if avgR(i) > 1.05*meanavgR | avgR(i) < 0.95*meanavgR
         avgR(i) = meanavgR;
     end
 end

 for i = 1:nfiles
     if avgG(i) > 1.05*meanavgG | avgG(i) < 0.95*meanavgG
        avgG(i) = meanavgG;
     end
 end

 for i = 1:nfiles
     if avgB(i) > 1.05*meanavgB | avgB(i) < 0.95*meanavgB
         avgB(i) = meanavgB;
     end
 end
%
xaxis = 1:10;
figure(1)
hold on
plot(xaxis,avgR,'r')
plot(xaxis,avgG,'g')
plot(xaxis,avgB,'b')
title('FPS = 70')
xlabel('# images')
ylabel('not normalized intensity')
%%
normavgR = 100.*(avgR./max(avgR));
normavgG = 100.*(avgG./max(avgG));
normavgB = 100.*(avgB./max(avgB));
%%
figure(2)
hold on
plot(xaxis,normavgR,'r')
plot(xaxis,normavgG,'g')
plot(xaxis,normavgB,'b')
ylim([90 100])
title('FPS = 70')
%% 
xlabel('# images')
ylabel('normalized intensity')
saveas(figure(2),[pwd '/FPS_IMAGES/70FPS.fig']);
%%
%%
% for ii=1:nfiles
%    current_filename = imagefiles(ii).name;
%    [Y, M, D, H, MN, S] = datevec(imagefiles(ii).datenum);
%    current_image = imread(current_filename);
%    current_image_crop = current_image(x_crop_0nm:x_crop_0nm+63,y_crop_0nm:y_crop_0nm+63);
%    current_image_avg_intensity = mean(mean(current_image_crop));
%    intensities(1,ii) = current_image_avg_intensity;
%    images{ii} = current_image_crop;
%    timestamps{ii} = [H, MN, S];
% end

% tot_time = timestamps{nfiles}(1,:) - timestamps{1}(1,:); %total time between first and last image
% t_in_secs = (tot_time(1,1)*3600) + (tot_time(1,2)*60) + tot_time(1,3);
% 
% intensities = 100.*(intensities./max(intensities));

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

%%
% x = linspace(0, t_in_secs, nfiles);
% y = intensities;
% plot(x,y)
% ylim([90 100])
% title('Intensity vs Time for Red Image')
toc