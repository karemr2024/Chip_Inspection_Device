function [maskedImage] = AreaSelection_Circle_Mod(originalImage,x,y)

% Demo to have the user click and draw a circle over an image, then blacken outside the circle and crop out the circular portion into a new image.
clc;% Clear the command window.
fprintf('Beginning to run %s.m ...\n', mfilename);
close all;  % Close all figures (except those of imtool.)
imtool close all;  % Close all imtool figures.
workspace;  % Make sure the workspace panel is showing
fontSize = 15;

% Get image.
[rows, columns, ~] = size(originalImage);

% Show circle over image.
figure(2)
imshow(originalImage);
axis('on', 'image');
hold on;
plot(x, y, 'r-', 'LineWidth', 2);
title('Original image with circle mask overlaid', 'FontSize', fontSize);
% Get a mask of the circle
mask = poly2mask(x, y, rows, columns);
% figure(3)
% imshow(mask);
% axis('on', 'image');
% title('Circle Mask', 'FontSize', fontSize);

% Mask the image with the circle.
% if numberOfColorChannels == 1
% 	maskedImage = originalImage; % Initialize with the entire image.
% 	maskedImage(~circleImage) = 0; % Zero image outside the circle mask.
% else
	% Mask the image.
	maskedImage = bsxfun(@times, originalImage, cast(mask, class(originalImage)));
% end

% Crop the image to the bounding box.
props = regionprops(mask, 'BoundingBox');
maskedImage = imcrop(maskedImage, props.BoundingBox);

% Display it in the lower right plot.
figure(3);
imshow(maskedImage, []);
% Change imshow to image() if you don't have the Image Processing Toolbox.
title('Image masked with the circle', 'FontSize', fontSize);
fprintf('Done running %s.m ...\n', mfilename);

end