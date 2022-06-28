tic

clc; clearvars; close all;

load("Imaging_Data.mat")
Ti_nData = table2array(readtable("Palm.csv"));
Ti_nData = Ti_nData(110:286,:);
Ti_n = interp(Ti_nData(:,2),15)'; %refractive index values of Ti
%
L_Ti = linspace(0.05,0.1,2655); %Silicon Oxide thickness on Ti from 50 to 100 nm
lambda_Ti = interp(Ti_nData(:,1),15)'; %lambda values from 400 to 680 nm

nsqrSiO =sellmeier(B_SiO,C_SiO,lambda_Ti);
nrSiO = (sqrt(nsqrSiO)+conj(sqrt(nsqrSiO)))/2;
%
Z1_Ti = [];
Gamma1_Ti = [];
%
for i = 1:numel(lambda_Ti)
for j = 1:numel(L_Ti)
[Gamma1_Ti(i,j),Z1_Ti(i,j)] = multidiels([1;nrSiO(1,i);Ti_n(1,i)],L_Ti(j).*nrSiO(1,i),lambda_Ti(1,i));
end
end
%
[MLam_Ti, MThicc_Ti] = meshgrid(lambda_Ti,L_Ti);
Gamma_Ti = conj(Gamma1_Ti).*Gamma1_Ti; %Multiply Gamma with conjugate to get rid of imaginary component

%%
figure(1)

surgraph = surf(MLam_Ti,MThicc_Ti,Gamma_Ti,'EdgeColor','none');
title('Ti/SiO2 Reflectance')
xlim([0.4 0.68])
xticks([0.4 0.44 0.48 0.52 0.56 0.60 0.64 0.68])
xticklabels([0.05 0.0571 0.0643 0.0714 0.0785 0.0856 0.0927 0.1])
ylim([0.05 0.1])
yticks([0.05 0.0571 0.0643 0.0714 0.0785 0.0856 0.0927 0.1])
yticklabels([0.4 0.44 0.48 0.52 0.56 0.60 0.64 0.68])
cb = colorbar;
cb.Location = 'eastoutside';
xlabel('L (\mum)');
ylabel('Lambda (\mum)');
zlabel('Reflectance (%)');
caxis([min(min(Gamma_Ti)) max(max(Gamma_Ti))]);

%

figure(2) % Thickness vs Reflectance for RGB wavelengths

hold on
plot(L_Ti,Gamma_Ti(valtoindex_lambda_Ti(cw_b),:),'b','LineWidth',2) %Light wavelength in blue
plot(L_Ti,Gamma_Ti(valtoindex_lambda_Ti(cw_g),:),'g','LineWidth',2) %Light wavelength in green
plot(L_Ti,Gamma_Ti(valtoindex_lambda_Ti(cw_r),:),'r','LineWidth',2) %Light wavelength in red
xlabel('Thickness (\mum)')
ylabel('Reflectance (%)')
title('Thickness vs Reflectance for RGB wavelengths')
legend('Blue', 'Green', 'Red')

figure(3) %Wavelength vs Reflectance for different thicknesses

hold on
plot(lambda_Ti,Gamma_Ti(:,valtoindex_L_Ti(0.05)),'r','LineWidth',2) %SiO2 Thickness at 50 nm 
plot(lambda_Ti,Gamma_Ti(:,valtoindex_L_Ti(0.1)),'g','LineWidth',2) %SiO2 Thickness at 100 nm

xlim([0.4 0.68])
xlabel('Wavelength \mum)')
ylabel('Reflectance (%)')
title('Wavelength vs Reflectance for 50 and 100 nm')
legend('L = 50 nm', 'L = 100 nm')

% Gamma matrix is proven. Gamma curves at lambda values 449, 521 nm AND L values 72,77,82 nm are 
% very similar to the curves seen in source graphs.

%%

%Now goal is to find reflectivity curves for red, green and blue. Center
%wavelength values are taken from the center wavelength values of camera

figure(4)

hold on
plot(L_Ti,Gamma_Ti(valtoindex_lambda_Ti(cw_b),:),'b','LineWidth',2) %Reflectivity curve at 460 nm (blue)
plot(L_Ti,Gamma_Ti(valtoindex_lambda_Ti(cw_g),:),'g','LineWidth',2) %Reflectivity curve at 530 nm (green)
plot(L_Ti,Gamma_Ti(valtoindex_lambda_Ti(cw_r),:),'r','LineWidth',2) %Reflectivity curve at 625 nm (red)
legend('Blue','Green','Red')
title('RGB Light Reflectivity')
xlabel('L (\mum)','FontSize',16);
ylabel('Reflectivity (%)','FontSize',16');
xlim([0.05 0.1])
hold off

%Spectrum Data for red, green and blue LEDs are shown in Figure 5.
%%
load("spectra.mat",'spectra') %Load spectrum data
blue = spectra{1,1}; %LED Spectrum for B
green = spectra{1,2}; %LED Spectrum for G
red = spectra{1,4}; %LED Spectrum for R
bluespectrum = blue(:,2)*100;
bluespectrum = bluespectrum(251:1620); %Blue Spectrum from 400 to 680 nm
bluespectrum = interp(bluespectrum,2);
greenspectrum = green(:,2)*100;
greenspectrum = greenspectrum(251:1620);
greenspectrum = interp(greenspectrum,2);%Green Spectrum from 400 to 680 nm
redspectrum = red(:,2)*100;
redspectrum = redspectrum(251:1620);
redspectrum = interp(redspectrum,2);%Red Spectrum from 400 to 680 nm

figure(5) %Plot RGB Spectrum

hold on
plot(lambda_Ti,bluespectrum,'b','LineWidth',2)
plot(lambda_Ti,greenspectrum,'g','LineWidth',2)
plot(lambda_Ti,redspectrum,'r','LineWidth',2)
title('LED Light Intensities')
xlim([0.4 0.68])
xlabel('Wavelength (nm)')
ylabel('Intensity')
hold off
%%
%Multiply LED spectrum with reflectivity curves to find reflected Intensity
for i = 1:2740
Ref_spec_red = Gamma(i,valtoindex_L(0))'.*redspectrum; %Reflectivity spectrum for red at 0 nm (R_r) 
Ref_spec_green = Gamma(i,valtoindex_L(0))'.*greenspectrum; %Reflectivity spectrum for green at 0 nm (R_g) 
Ref_spec_blue = Gamma(i,valtoindex_L(0))'.*bluespectrum; %Reflectivity spectrum for blue at 0 nm (R_b)
end 

reflectivity_curve_0nm = Gamma(:,valtoindex_L(0)); %Reflectivity curve at 0 nm (used above)

%Plot RGB reflectivity curves for Si chip with no oxide (SiO2 thickness = 0 nm)

figure(6)
hold on
plot(lambda,reflectivity_curve_0nm,'k-','LineWidth',2) %Reflectivity curve at 0 nm 
plot(lambda,redspectrum,'r','Linewidth',2) %Red spectrum before multiplication
plot(lambda,greenspectrum,'g','Linewidth',2) %Green spectrum before multiplication
plot(lambda,bluespectrum,'b','Linewidth',2) %Blue spectrum before multiplication
plot(lambda,Ref_spec_red,'r','Linewidth',2) %Red spectrum after multiplication
plot(lambda,Ref_spec_green,'g','Linewidth',2) %Green spectrum after multiplication
plot(lambda,Ref_spec_blue,'b','Linewidth',2) %Blue spectrum after multiplication
xlabel('Wavelength (nm)')
xlim([0.4 0.68])
ylabel('Reflectivity')
title('Reflectivity Curve for Silicon at 0 nm (no oxide)')
legend('Reflectivity Curve for Silicon at 0 nm','Red Spectrum','Green Spectrum','Blue Spectrum','location','north')

%Find FWHM of new multiplied images
halfMaxb = (min(Ref_spec_blue) + max(Ref_spec_blue)) / 2;
index1b = find(Ref_spec_blue >= halfMaxb, 1, 'first');
index2b = find(Ref_spec_blue >= halfMaxb, 1, 'last');
fwhmb = index2b-index1b + 1;
fwhmxb = lambda(index2b) - lambda(index1b);

halfMaxg = (min(Ref_spec_green) + max(Ref_spec_green)) / 2;
index1g = find(Ref_spec_green >= halfMaxg, 1, 'first');
index2g = find(Ref_spec_green >= halfMaxg, 1, 'last');
fwhmg = index2g-index1g + 1;
fwhmxg = lambda(index2g) - lambda(index1g);

halfMaxr = (min(Ref_spec_red) + max(Ref_spec_red)) / 2;
index1r = find(Ref_spec_red >= halfMaxr, 1, 'first');
index2r = find(Ref_spec_red >= halfMaxr, 1, 'last');
fwhmr = index2r-index1r + 1;
fwhmxr = lambda(index2r) - lambda(index1r);

%Find center wavelength of I_Out
cw_b = (lambda(index2b)+lambda(index1b))/2;
cw_g = (lambda(index2g)+lambda(index1g))/2;
cw_r = (lambda(index2r)+lambda(index1r))/2;

%Find reflectance values at center wavelength
Ref_at_blue = reflectivity_curve_0nm(valtoindex_lambda(cw_b)); %Reflectivity at blue (455 nm)
Ref_at_green = reflectivity_curve_0nm(valtoindex_lambda(cw_g)); %Reflectivity at green (521.2 nm)
Ref_at_red = reflectivity_curve_0nm(valtoindex_lambda(cw_r)); %Reflectivity at red (634.2 nm)



% %Find intersection points of reflectivity curve and RGB spectrum
% 
% intsct_b1 = find(bluespectrum >= reflectivity_curve_0nm, 1, 'first');
% intsct_b2 = find(bluespectrum >= reflectivity_curve_0nm, 1, 'last');
% intsct_g1 = find(greenspectrum >= reflectivity_curve_0nm, 1, 'first')
% intsct_g2 = find(greenspectrum >= reflectivity_curve_0nm, 1, 'last')
% intsct_r1 = find(redspectrum >= reflectivity_curve_0nm, 1, 'first')
% intsct_r2 = find(redspectrum >= reflectivity_curve_0nm, 1, 'last')

figure(7)

I_refred = []; %Reflected Intensity for red
I_refgreen = []; %Reflected Intensity for green
I_refblue = []; %Reflected Intensity for blue

for i=1:2740
    I_refred(:,i) = sum(Gamma(:,i).*Ref_spec_red);          %Multiply reflectance with I_Out for every thickness and sum to find total reflected intensity
    I_refgreen(:,i) = sum(Gamma(:,i).*Ref_spec_green);      %Multiply reflectance with I_Out for every thickness and sum to find total reflected intensity
    I_refblue(:,i) = sum(Gamma(:,i).*Ref_spec_blue);        %Multiply reflectance with I_Out for every thickness and sum to find total reflected intensity
end 

hold on
plot(L,I_refblue/100,'b','LineWidth',2)
plot(L,I_refgreen/100,'g','LineWidth',2)
plot(L,I_refred/100,'r','LineWidth',2)
title('Reflected Intensity (I_O_u_t) for RGB from 0 nm to 140 nm')
legend('Blue','Green','Red')
xlabel('L (\mum)','FontSize',16);
ylabel('Reflected Intensity (I_O_u_t)','FontSize',16);
hold off

figure(10) %Simulated color for 0-300 nm
I_tot = cat(3, I_refred, I_refgreen, I_refblue);
I_totnorm = (I_tot./max(I_tot)).*255;

for i = 1:6
    I_totnorm = [I_totnorm; I_totnorm];
end

imtot = uint8(round(I_totnorm));
imagesc(imtot)
set(gca,'YDir','normal')
ylim([1,30])
title('Simulated Color')
xlabel('Thickness')
xticks([valtoindex_L(0) valtoindex_L(0.05) valtoindex_L(0.10) valtoindex_L(0.15) valtoindex_L(0.20) valtoindex_L(0.25) valtoindex_L(0.30)])
xticklabels([0 0.05 0.10 0.15 0.20 0.25 0.30])

figure(11) %Color for 120 nm

im_120 = imtot(:,valtoindex_L(0.12),:);

for i = 1:6
    im_120 = [im_120 im_120];
end

imagesc(im_120)

toc 