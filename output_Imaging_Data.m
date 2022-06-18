clc; clearvars; close all;

% Define variables for input in multidiels function: 

L = linspace(0,0.3,2740); %SiO2 Thickness from 0 nm (Si) to 140 nm. 
lambda = linspace(0.4,0.68,2740); %Wavelength of LED from 400 nm to 680 nm. 
% Functions valtoindex_L & valtoindex_lambda are used to convert L or
% lambda value to MATLAB index

theta = 0; %incidence angle from left medium (in degrees)

%Below are B and C constants for Si and SiO2 to use in Sellmeier equation.

B_Si = [10.6684293 0.0030434748 1.54133408]; %B Constants for Si [B1, B2, B3]
B_SiO = [0.6961663 0.4079426 0.8974794]; %B Constants for SiO2 [B1, B2, B3]

C_Si = [0.301516485 1.13475115 1104]; %C Constants for Si [C1, C2, C3]
C_SiO = [0.0684043 0.1162414 9.896161]; %C Constants for SiO2 [C1, C2, C3]

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

for i = 1:2740
Ref_spec_red = Gamma(i,valtoindex_L(0))'.*redspectrum; %Reflectivity spectrum for red at 0 nm (R_r) 
Ref_spec_green = Gamma(i,valtoindex_L(0))'.*greenspectrum; %Reflectivity spectrum for green at 0 nm (R_g) 
Ref_spec_blue = Gamma(i,valtoindex_L(0))'.*bluespectrum; %Reflectivity spectrum for blue at 0 nm (R_b)
end 

reflectivity_curve_0nm = Gamma(:,valtoindex_L(0));

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

Ref_at_0nm_R = Gamma(valtoindex_lambda(cw_r),valtoindex_L(0));
Ref_at_0nm_G = Gamma(valtoindex_lambda(cw_g),valtoindex_L(0));
Ref_at_0nm_B = Gamma(valtoindex_lambda(cw_b),valtoindex_L(0));
Ref_at_120nm_R = Gamma(valtoindex_lambda(cw_r),valtoindex_L(0.12));
Ref_at_120nm_G = Gamma(valtoindex_lambda(cw_g),valtoindex_L(0.12));
Ref_at_120nm_B = Gamma(valtoindex_lambda(cw_b),valtoindex_L(0.12));

clearvars -except Gamma L lambda cw_b cw_g cw_r Ref_at_0nm_R Ref_at_0nm_G Ref_at_0nm_B Ref_at_120nm_R Ref_at_120nm_G Ref_at_120nm_B
save("Imaging_Data.mat")
