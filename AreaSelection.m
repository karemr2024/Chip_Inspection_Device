% function [I2,dataX,dataYY] = AreaSelection(b)
 GH = figure; imshow(R_Tiff_0nm_Crop)%%%Shows the input image%%%%
  waitforbuttonpress  %%%%%%Press left click on the image%%%%
  ss = size(R_Tiff_0nm_Crop);
  point1 = get(gcf,'CurrentPoint') 
  rect = [point1(1,1) point1(1,2) 64 64]; %Coordinates of rectangle%
  [r2] = dragrect(rect);                    %%%Drag the rectangle while keeping leftclick pressed and leave the click when region to be selected is decided%%
  [dataX, dataY] = pix2data(r2(1,1),r2(1,2));
  ggr = ss(1,1) - dataY ;
  dataY = ggr+0.5;
  dataX=round(dataX+0.5);                      %Top left hand side X coordinate of the image%
  dataYY =round(dataY-r2(1,4));                %Top left hand side Y coordinate of the image%  
  I2 = imcrop(R_Tiff_0nm_Crop,[dataX dataYY 63 63]);%Final Cropped image%
  close(GH)