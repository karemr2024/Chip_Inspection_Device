clc; clearvars; close all;

% Define variables for input in multidiels function: 

L = linspace(0,0.14,1370); %SiO2 Thickness from 0 nm (Si) to 140 nm. 
lambda = linspace(0.4,0.68,1370); %Wavelength of LED from 400 nm to 680 nm. 
% Functions valtoindex_L & valtoindex_lambda are used to convert L or
% lambda value to MATLAB index

theta = 0; %incidence angle from left medium (in degrees)

%Below are B and C constants for Si and SiO2 to use in Sellmeier equation.

B_Si = [10.6684293, 0.0030434748, 1.54133408]; %B Constants for Si [B1, B2, B3]
B_SiO = [0.6961663, 0.4079426, 0.8974794]; %B Constants for SiO2 [B1, B2, B3]

C_Si = [0.301516485, 1.13475115, 1104]; %C Constants for Si [C1, C2, C3]
C_SiO = [0.0684043, 0.1162414, 9.896161]; %C Constants for SiO2 [C1, C2, C3]

nsqrSi =sellmeier(B_Si,C_Si,lambda); %Find refractive index values for Si at each wavelength
nsqrSiO =sellmeier(B_SiO,C_SiO,lambda); %Find refractive index values for SiO2 at each wavelength

nrSi = (sqrt(nsqrSi)+conj(sqrt(nsqrSi)))/2; %Si refractive index to input in multidiels
nrSiO = (sqrt(nsqrSiO)+conj(sqrt(nsqrSiO)))/2; %SiO2 refractive index to input in multidiels

% Calculate and display Gamma matrix from inputs above. Gamma(lambda,L) is a matrix for Si/SiO2
% reflectivity. The variables are light wavelength (lambda) and SiO2 thickness (L). 
% Usage: [Gamma,Z] = multidiels(n,L,lambda,theta,pol) Theta and pol are 0, so they are not included

Z1 = [];
Gamma1 = [];

for i = 1:numel(lambda)
for j = 1:numel(L)
[Gamma1(i,j),Z1(i,j)] = multidiels([1;nrSiO(1,i);nrSi(1,i)],L(j).*nrSiO(1,i),lambda(1,i));
end
end

[MLam, MThicc] = meshgrid(lambda,L);
Gamma = conj(Gamma1).*Gamma1; %Multiply Gamma with conjugate to get rid of imaginary component
surgraph = surf(MLam,MThicc,Gamma,'EdgeColor','none');
title('Si/SiO2 Reflectance')
xlim([0.4 0.68])
ylim([0 0.14])
cb = colorbar;

cb.Location = 'eastoutside';
xlabel('Lambda (\mum)');
ylabel('L (\mum)');
zlabel('Reflectance (%)');
caxis([0 0.5]);

% PROOF
% Reflectivity from Gamma matrix at wavelengths 449 nm and 521 nm & thicknesses 72 nm, 77 nm, and 
% 82 nm are compared to figures 3 & 4 in https://www.researchgate.net/figure/Surface-reflectivity-versus-silicon-dioxide-layer-thickness_fig3_230952570

figure(2) % Thickness vs Reflectance for different wavelengths

lambda_449nm = valtoindex_lambda(0.449); %wavelength at 449 nm
lambda_521nm = valtoindex_lambda(0.521); %wavelength at 521 nm

hold on
plot(L,Gamma(lambda_449nm,:),'b','LineWidth',2)
plot(L,Gamma(lambda_521nm,:),'g','LineWidth',2)
xlabel('Thickness (\mum)')
ylabel('Reflectance (%)')
title('Thickness vs Reflectance for different wavelengths')
legend('449 nm', '521 nm')

figure(3) %Wavelength vs Reflectance for different thicknesses

L_72nm = valtoindex_L(0.072); %SiO2 Thickness at 72 nm 
L_77nm = valtoindex_L(0.077); %SiO2 Thickness at 77 nm 
L_82nm = valtoindex_L(0.082); %SiO2 Thickness at 82 nm 

hold on
plot(lambda,Gamma(:,L_72nm),'r','LineWidth',2)
plot(lambda,Gamma(:,L_77nm),'g','LineWidth',2)
plot(lambda,Gamma(:,L_82nm),'b','LineWidth',2)
xlabel('Wavelength \mum)')
ylabel('Reflectance (%)')
title('Wavelength vs Reflectance for different thicknesses')
legend('72 nm', '77 nm', '82 nm')

% Gamma matrix is proven. Gamma curves at lambda values 449, 521 nm AND L values 72,77,82 nm are 
% very similar to the curves seen in source graphs.

%Now goal is to find reflectivity curves for red, green and blue. Center
%wavelength values are taken from the center wavelength values of 

figure(4)

hold on
plot(L,Gamma(valtoindex_lambda(0.46),:),'b','LineWidth',2) %Reflectivity curve at 460 nm (blue)
plot(L,Gamma(valtoindex_lambda(0.53),:),'g','LineWidth',2) %Reflectivity curve at 530 nm (green)
plot(L,Gamma(valtoindex_lambda(0.625),:),'r','LineWidth',2) %Reflectivity curve at 625 nm (red)
legend('Blue','Green','Red')
title('RGB Light Reflectivity')
xlabel('L (\mum)','FontSize',16);
ylabel('Reflectivity (%)','FontSize',16');
xlim([0 0.14])
hold off

%we have to find intensity_in from here - figure out tomorrow

%Spectrum Data for red, green and blue LEDs are shown in Figure 5.

load("spectra.mat",'spectra') %Load spectrum data
blue = spectra{1,1}; %LED Spectrum for B
green = spectra{1,2}; %LED Spectrum for G
red = spectra{1,4}; %LED Spectrum for R
bluespectrum = blue(:,2)*100;
bluespectrum = bluespectrum(251:1620); %Blue Spectrum from 400 to 680 nm
greenspectrum = green(:,2)*100;
greenspectrum = greenspectrum(251:1620); %Green Spectrum from 400 to 680 nm
redspectrum = red(:,2)*100;
redspectrum = redspectrum(251:1620); %Red Spectrum from 400 to 680 nm

figure(5) %Plot RGB Spectrum

hold on
plot(lambda,bluespectrum,'b','LineWidth',2)
plot(lambda,greenspectrum,'g','LineWidth',2)
plot(lambda,redspectrum,'r','LineWidth',2)
title('LED Light Intensities')
xlabel('Wavelength (nm)')
ylabel('Reflectivity')
hold off

%Multiply LED spectrum with reflectivity curves to find to actual spectrum data 
for i = 1:1370
real_redspectrum = Gamma(i,:)'.*redspectrum; 
real_greenspectrum = Gamma(i,:)'.*greenspectrum; 
real_bluespectrum = Gamma(i,:)'.*bluespectrum; 
end 

real_redspectrum = Gamma(valtoindex_L(0),:)'.*redspectrum; %More accurate spectrum data for R LED
real_greenspectrum = Gamma(valtoindex_L(0),:)'.*greenspectrum; %More accurate spectrum data for G LED
real_bluespectrum = Gamma(valtoindex_L(0),:)'.*bluespectrum; %More accurate spectrum data for B LED

%Show RGB reflectivity curves for Si chip with no oxide (SiO2 thickness = 0 nm)

figure(6)
hold on
plot(lambda,Gamma(valtoindex_L(0),:),'k-','LineWidth',2)
plot(lambda,redspectrum,'r','Linewidth',2)
plot(lambda,greenspectrum,'g','Linewidth',2)
plot(lambda,bluespectrum,'b','Linewidth',2)
plot(lambda,real_redspectrum,'r','Linewidth',2)
plot(lambda,real_greenspectrum,'g','Linewidth',2)
plot(lambda,real_bluespectrum,'b','Linewidth',2)
xlabel('Wavelength (nm)')
ylabel('Reflectivity')
title('RGB Reflectivity Curves for Silicon (no oxide)')
legend('Reflectivity Curve','Red Spectrum','Green Spectrum','Blue Spectrum','location','north')

figure(7)










