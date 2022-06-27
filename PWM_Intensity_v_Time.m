tic
% Get list of all BMP files in this directory
% DIR returns as a structure array.  You will need to use () and . to get
% the file names.
% Only works for files taken within the same day
clearvars; clc; close all;
%START NEW CODE

tiff_info_R = imfinfo('StackR.tiff'); % return tiff structure, one element per image
tiff_stack_R = imread('StackR.tiff', 1) ; % read in first image
tiff_info_G = imfinfo('StackG.tiff'); % return tiff structure, one element per image
tiff_stack_G = imread('StackG.tiff', 1) ; % read in first image
tiff_info_B = imfinfo('StackB.tiff'); % return tiff structure, one element per image
tiff_stack_B = imread('StackB.tiff', 1) ; % read in first image
%concatenate each successive tiff to tiff_stack
for i = 2 : size(tiff_info_R, 1)
    temp_tiff_R = imread('StackR.tiff', i);
    tiff_stack_R = cat(3 , tiff_stack_R, temp_tiff_R);
    temp_tiff_G = imread('StackG.tiff', i);
    tiff_stack_G = cat(3 , tiff_stack_G, temp_tiff_G);
    temp_tiff_B = imread('StackB.tiff', i);
    tiff_stack_B = cat(3 , tiff_stack_B, temp_tiff_B);
end

tiff_stack_R = double(tiff_stack_R);
tiff_stack_G = double(tiff_stack_G);
tiff_stack_B = double(tiff_stack_B);
%%

%END NEW CODE
% imagefiles = dir('**//*.tif'); 
%%
size = size(tiff_stack_B);
nfiles = size(3);
%%

% intensities = zeros(1,nfiles);

% Sort images according to date and time
% imagefiles = imagefiles(~[imagefiles.isdir]);
% [~,idx] = sort([imagefiles.datenum]);
% imagefiles = imagefiles(idx);

% [display_image_crop,x_crop_0nm,y_crop_0nm] = AreaSelection(imread(imagefiles(2).name));

for i = 1:nfiles
    avgR(i) = mean(mean(tiff_stack_R(:,:,i)));   
    avgG(i) = mean(mean(tiff_stack_G(:,:,i)));   
    avgB(i) = mean(mean(tiff_stack_B(:,:,i)));   
end

%%
normavgR = 100.*(avgR./max(avgR));
normavgG = 100.*(avgG./max(avgG));
normavgB = 100.*(avgB./max(avgB));
%%
time = 20
xaxis = linspace(1,time,nfiles)
figure(1)
hold on
plot(xaxis,normavgR,'r')
plot(xaxis,normavgG,'g')
plot(xaxis,normavgB,'b')
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