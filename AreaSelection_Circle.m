function [maskedImage,x,y,point1] = AreaSelection_Circle(originalImage)

% Demo to have the user click and draw a circle over an image, then blacken outside the circle and crop out the circular portion into a new image.
clc;% Clear the command window.
fprintf('Beginning to run %s.m ...\n', mfilename);
close all;  % Close all figures (except those of imtool.)
imtool close all;  % Close all imtool figures.
workspace;  % Make sure the workspace panel is showing
fontSize = 15;
% Get image.
[rows, columns, ~] = size(originalImage);
fig = figure
imshow(originalImage);
axis('on', 'image');
% Maximize the window to make it easier to draw.
g = gcf;
g.WindowState = 'maximized';

uiwait(helpdlg('Select ROI enclosing the spot.'));

h.Radius = 0;
while h.Radius == 0
	h = drawcircle('Color','k','FaceAlpha',0);
	if h.Radius == 0
		uiwait(helpdlg('You double-clicked.  You need to single click, then drag, then single click again.'));
	end
end
% Get coordinates of the circle.,'DrawingArea',[x,y,w,h]
angles = linspace(0, 2*pi, 10000);
x = cos(angles) * h.Radius + h.Center(1);
y = sin(angles) * h.Radius + h.Center(2);
% Get a mask of the circle
mask = poly2mask(x, y, rows, columns);

% Mask the image.
	maskedImage = bsxfun(@times, originalImage, cast(mask, class(originalImage)));
% end

% Crop the image to the bounding box.
props = regionprops(mask, 'BoundingBox');
maskedImage = imcrop(maskedImage, props.BoundingBox);

fprintf('Done running %s.m ...\n', mfilename);

end