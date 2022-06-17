% Get list of all BMP files in this directory
% DIR returns as a structure array.  You will need to use () and . to get
% the file names.
% Only works for files taken within the same day

imagefiles = dir('**//*.tif');      
nfiles = length(imagefiles);
intensities = zeros(1,nfiles);

% Sort images according to date and time
imagefiles = imagefiles(~[imagefiles.isdir]);
[~,idx] = sort([imagefiles.datenum]);
imagefiles = imagefiles(idx);

[display_image_crop,x_crop_0nm,y_crop_0nm] = AreaSelection(imread(imagefiles(1).name));

for ii=1:nfiles
   current_filename = imagefiles(ii).name;
   [Y, M, D, H, MN, S] = datevec(imagefiles(ii).datenum);
   current_image = imread(current_filename);
   current_image_crop = current_image(x_crop_0nm:x_crop_0nm+63,y_crop_0nm:y_crop_0nm+63);
   current_image_avg_intensity = mean(mean(current_image_crop));
   intensities(1,ii) = current_image_avg_intensity;
   images{ii} = current_image_crop;
   timestamps{ii} = [H, MN, S];
end

tot_time = timestamps{nfiles}(1,:) - timestamps{1}(1,:); %total time between first and last image
t_in_secs = (tot_time(1,1)*3600) + (tot_time(1,2)*60) + tot_time(1,3);

x = linspace(0, t_in_secs, nfiles);
y = intensities;
plot(x,y)
