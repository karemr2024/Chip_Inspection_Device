%% START

tic
clc; clearvars; close all;
fprintf('Beginning to run %s.m ...\n', mfilename);

%% Define variables for input in multidiels function:

numval = 2740;
L = linspace(0,0.3,numval); %SiO2 Thickness from 0 nm (Si) to 300 nm. 
lambda = linspace(0.4,0.68,numval); %Wavelength of LED from 400 nm to 680 nm. 

%% Below are B and C constants for Si and SiO2 to use in Sellmeier equation. 
%% B and C values are taken from https://refractiveindex.info/

B_Si = [10.6684293 0.0030434748 1.54133408]; %B Constants for Si [B1, B2, B3]
B_SiO = [0.6961663 0.4079426 0.8974794]; %B Constants for SiO2 [B1, B2, B3]

C_Si = [0.301516485 1.13475115 1104]; %C Constants for Si [C1, C2, C3]
C_SiO = [0.0684043 0.1162414 9.896161]; %C Constants for SiO2 [C1, C2, C3]

%% Calculate refractive index values for Si and SiO2

nsqrSi =sellmeier(B_Si,C_Si,lambda); %Find refractive index values for Si at each wavelength
nsqrSiO =sellmeier(B_SiO,C_SiO,lambda); %Find refractive index values for SiO2 at each wavelength

nrSi = (sqrt(nsqrSi)+conj(sqrt(nsqrSi)))/2; %Si refractive index to input in multidiels
nrSiO = (sqrt(nsqrSiO)+conj(sqrt(nsqrSiO)))/2; %SiO2 refractive index to input in multidiels



%% Calculate Gamma matrix from inputs above. Gamma(lambda,L) is a matrix for Si/SiO2
%% reflectivity. The variables are light wavelength (lambda) and SiO2 thickness (L). 
%% Usage: [Gamma,Z] = multidiels(n,L,lambda,theta,pol) Theta and pol are 0, so they are not included.

Z1 = [];
Gamma1 = [];

for i = 1:numval
for j = 1:numval
[Gamma1(i,j),Z1(i,j)] = multidiels([1;nrSiO(1,i);nrSi(1,i)],L(j).*nrSiO(1,i),lambda(1,i));
end
end

Gamma = conj(Gamma1).*Gamma1; %Multiply Gamma with conjugate to get rid of imaginary component

%% Display Gamma matrix

figure(1)

surf(L,lambda,Gamma,'EdgeColor','none');
title('Si/SiO2 Reflectance')
xlim([0 0.3])
ylim([0.4 0.68])
xlabel('L (\mum)');
ylabel('Lambda (\mum)');
zlabel('Reflectance (%)');
cb = colorbar;
cb.Location = 'eastoutside';
caxis([0 0.5])

%% Find reflectivity curves for red, green and blue. Wavelength values
%% are taken from the center wavelength values of camera

figure(2)

hold on
plot(L,Gamma(valtoindex_lambda(0.46),:),'b','LineWidth',2) %Reflectivity curve at 460 nm (blue)
plot(L,Gamma(valtoindex_lambda(0.53),:),'g','LineWidth',2) %Reflectivity curve at 530 nm (green)
plot(L,Gamma(valtoindex_lambda(0.625),:),'r','LineWidth',2) %Reflectivity curve at 625 nm (red)
legend('Blue','Green','Red')
title('RGB Light Reflectivity')
xlabel('L (\mum)','FontSize',16);
ylabel('Reflectivity (%)','FontSize',16');
xlim([0 0.3])

%% Spectrum Data for red, green and blue LEDs are shown in Figure 5.

load("spectra.mat",'spectra') %Load spectrum data
spectra = cell2mat(spectra);
spectra = spectra(251:1620,:);
bluespectrum = interp(spectra(:,2)*100,2);  %Blue Spectrum from 400 to 680 nm
greenspectrum = interp(spectra(:,4)*100,2); %Green Spectrum from 400 to 680 nm 
redspectrum = interp(spectra(:,8)*100,2);   %Red Spectrum from 400 to 680 nm 

bluespectrum = bluespectrum./max(max(bluespectrum));
greenspectrum = greenspectrum./max(max(greenspectrum));
redspectrum = redspectrum./max(max(redspectrum)); 

%% Plot RGB Spectrum

figure(3)

hold on
plot(lambda,bluespectrum,'b','LineWidth',2)
plot(lambda,greenspectrum,'g','LineWidth',2)
plot(lambda,redspectrum,'r','LineWidth',2)
title('LED Light Intensities')
xlim([0.4 0.68])
xlabel('Wavelength (nm)')
ylabel('Intensity')

%% Multiply LED spectrum with reflectivity curves to find reflected Intensity

reflectivity_curve_0nm = Gamma(:,valtoindex_L(0)); %Reflectivity curve at 0 nm (used below)

for i = 1:numval
Ref_spec_red = Gamma(i,valtoindex_L(0))'.*redspectrum; %Reflectivity spectrum for red at 0 nm (R_r) 
Ref_spec_green = Gamma(i,valtoindex_L(0))'.*greenspectrum; %Reflectivity spectrum for green at 0 nm (R_g) 
Ref_spec_blue = Gamma(i,valtoindex_L(0))'.*bluespectrum; %Reflectivity spectrum for blue at 0 nm (R_b)
end 

%% Plot RGB reflectivity curves for Si chip with no oxide (SiO2 thickness = 0 nm)

figure(4)
hold on

plot(lambda,reflectivity_curve_0nm,'k-','LineWidth',2) %Reflectivity curve at 0 nm 

plot(lambda,redspectrum,'r','Linewidth',2) %Red spectrum before multiplication
plot(lambda,greenspectrum,'g','Linewidth',2) %Green spectrum before multiplication
plot(lambda,bluespectrum,'b','Linewidth',2) %Blue spectrum before multiplication

plot(lambda,Ref_spec_red,'r','Linewidth',2) %Red spectrum after multiplication
plot(lambda,Ref_spec_green,'g','Linewidth',2) %Green spectrum after multiplication
plot(lambda,Ref_spec_blue,'b','Linewidth',2) %Blue spectrum after multiplication

xlim([0.4 0.68])
xlabel('Wavelength (nm)')
ylabel('Reflectivity')
title('Reflectivity Curve for Silicon at 0 nm (no oxide)')
legend('Reflectivity Curve for Silicon at 0 nm','Red Spectrum and I_O_u_t','Green Spectrum and I_O_u_t','Blue Spectrum and I_O_u_t','location','north')

%% Find FWHM of new multiplied curves

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

%% Find center wavelength of I_Out

cw_b = (lambda(index2b)+lambda(index1b))/2;
cw_g = (lambda(index2g)+lambda(index1g))/2;
cw_r = (lambda(index2r)+lambda(index1r))/2;

%% Find reflectance values at center wavelength

Ref_at_blue = reflectivity_curve_0nm(valtoindex_lambda(cw_b)); %Reflectivity at blue (455 nm)
Ref_at_green = reflectivity_curve_0nm(valtoindex_lambda(cw_g)); %Reflectivity at green (521.2 nm)
Ref_at_red = reflectivity_curve_0nm(valtoindex_lambda(cw_r)); %Reflectivity at red (634.2 nm)

%% Calculate total reflected intensity by multiplying reflectance with reflected intensity (I_Out)

I_refred = []; %Reflected Intensity for red
I_refgreen = []; %Reflected Intensity for green
I_refblue = []; %Reflected Intensity for blue

for i=1:numval
    I_refred(:,i) = sum(Gamma(:,i).*Ref_spec_red);       
    I_refgreen(:,i) = sum(Gamma(:,i).*Ref_spec_green);      
    I_refblue(:,i) = sum(Gamma(:,i).*Ref_spec_blue);        
end 

%% Plot total reflected intensity for RGB

figure(5)

hold on
plot(L,I_refblue/100,'b','LineWidth',2)
plot(L,I_refgreen/100,'g','LineWidth',2)
plot(L,I_refred/100,'r','LineWidth',2)
title('Reflected Intensity (I_O_u_t) for RGB from 0 nm to 300 nm')
legend('Blue','Green','Red')
xlabel('L (\mum)','FontSize',16);
ylabel('Reflected Intensity (I_O_u_t)','FontSize',16);

%% Display simulated color for 0-300 nm by comparing RGB reflected intensity values

figure(6)

I_tot = cat(3, I_refred, I_refgreen, I_refblue);
I_totnorm = (I_tot./max(I_tot)).*255;

for i = 1:6                                
    I_totnorm = [I_totnorm; I_totnorm]; %Add to each other for visualization purposes
end

imtot = uint8(round(I_totnorm));
imagesc(imtot)
set(gca,'YDir','normal')
ylim([1,30])
title('Simulated Color Si-SiO_2')
xlabel('Thickness')
xticks([valtoindex_L(0) valtoindex_L(0.05) valtoindex_L(0.10) valtoindex_L(0.15) valtoindex_L(0.20) valtoindex_L(0.25) valtoindex_L(0.30)])
xticklabels([0 0.05 0.10 0.15 0.20 0.25 0.30])

%% Color for 120 nm

figure(7) 

im_120 = imtot(:,valtoindex_L(0.12),:);

for i = 1:6
    im_120 = [im_120 im_120];
end

imagesc(im_120)

%% Additional variables needed for data

Gamma_Si = abs((1 - nsqrSi)./ (1 + nsqrSi)).^2 ;

Gamma_B = Gamma/Gamma_Si(valtoindex_lambda(cw_b));
Gamma_G = Gamma/Gamma_Si(valtoindex_lambda(cw_g));
Gamma_R = Gamma/Gamma_Si(valtoindex_lambda(cw_r));

Ref_at_0nm_R = Gamma(valtoindex_lambda(cw_r),valtoindex_L(0));
Ref_at_0nm_G = Gamma(valtoindex_lambda(cw_g),valtoindex_L(0));
Ref_at_0nm_B = Gamma(valtoindex_lambda(cw_b),valtoindex_L(0));
Ref_at_120nm_R = Gamma(valtoindex_lambda(cw_r),valtoindex_L(0.12));
Ref_at_120nm_G = Gamma(valtoindex_lambda(cw_g),valtoindex_L(0.12));
Ref_at_120nm_B = Gamma(valtoindex_lambda(cw_b),valtoindex_L(0.12));

%% Save variables needed in other scripts

save Simulation_Data.mat

%% PROOF AND COMPARISON OF GAMMA MATRIX WITH SOURCE
%% 
%% Reflectivity from Gamma matrix at wavelengths 449 nm and 521 nm & thicknesses 72 nm, 77 nm, and 
%% 82 nm are compared to figures 3 & 4 in https://www.researchgate.net/figure/Surface-reflectivity-versus-silicon-dioxide-layer-thickness_fig3_230952570
%%
%% UN/COMMENT BELOW HERE

% figure % Thickness vs Reflectance for different wavelengths
% hold on
% plot(L,Gamma(valtoindex_lambda(0.449),:),'b','LineWidth',2) %Light wavelength at 449 nm
% plot(L,Gamma(valtoindex_lambda(0.521),:),'g','LineWidth',2) %Light wavelength at 521 nm
% xlabel('Thickness (\mum)')
% ylabel('Reflectance (%)')
% title('Thickness vs Reflectance for different wavelengths')
% legend('449 nm', '521 nm')
% 
% figure %Wavelength vs Reflectance for different thicknesses
% 
% hold on
% plot(lambda,Gamma(:,valtoindex_L(0.072)),'r','LineWidth',2) %SiO2 Thickness at 72 nm 
% plot(lambda,Gamma(:,valtoindex_L(0.077)),'g','LineWidth',2) %SiO2 Thickness at 77 nm
% plot(lambda,Gamma(:,valtoindex_L(0.082)),'b','LineWidth',2) %SiO2 Thickness at 82 nm
% xlim([0.4 0.68])
% xlabel('Wavelength \mum)')
% ylabel('Reflectance (%)')
% title('Wavelength vs Reflectance for different thicknesses')
% legend('L = 72 nm', 'L = 77 nm', 'L = 82 nm')
% 
% %Gamma matrix is proven. Gamma curves at lambda values 449, 521 nm AND L values 72,77,82 nm are 
% %very similar to the curves seen in source graphs.

%% UN/COMMENT ENDS

toc 