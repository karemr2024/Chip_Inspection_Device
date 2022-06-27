tic
   % The program section to time. 

clc; clearvars; close all;

load("Osram_Data.mat")
Gamma = Gamma(1:2691,1:2691);
L = L(1:2691);
lambda = lambda(1:2691);
Osram_Led_Spec_NoInt = table2array(readtable("LED_Engin_Spectrum_Data.csv"));

for i=1:23
Osram_Led_Spec(:,i) = interp(Osram_Led_Spec_NoInt(:,i),3); 
end

Osram_lambda = Osram_Led_Spec(:,1); %Wavelength from 350nm to 1100 nm

%Spectrum Datas (Normalized to 1)

Spec_Warm_White_2200K = Osram_Led_Spec(:,2)/(max(max(Osram_Led_Spec(:,2))));
Spec_Warm_White_3000K = Osram_Led_Spec(:,3)/(max(max(Osram_Led_Spec(:,3))));
Spec_Neutral_White_4000K = Osram_Led_Spec(:,4)/(max(max(Osram_Led_Spec(:,4))));
Spec_Cool_White_5500K = Osram_Led_Spec(:,5)/(max(max(Osram_Led_Spec(:,5))));
Spec_Cool_White_6500K = Osram_Led_Spec(:,6)/(max(max(Osram_Led_Spec(:,6))));
Spec_UVA_365nm = Osram_Led_Spec(:,7)/(max(max(Osram_Led_Spec(:,7))));
Spec_Violet_385nm = Osram_Led_Spec(:,8)/(max(max(Osram_Led_Spec(:,8))));
Spec_Violet_395nm = Osram_Led_Spec(:,9)/(max(max(Osram_Led_Spec(:,9))));
Spec_Violet_405nm = Osram_Led_Spec(:,10)/(max(max(Osram_Led_Spec(:,10))));
Spec_DeepBlue_436nm = Osram_Led_Spec(:,11)/(max(max(Osram_Led_Spec(:,11))));
Spec_RoyalBlue_453nm = Osram_Led_Spec(:,12)/(max(max(Osram_Led_Spec(:,12))));
Spec_DentalBlue_460nm = Osram_Led_Spec(:,13)/(max(max(Osram_Led_Spec(:,13))));
Spec_Cyan_500nm = Osram_Led_Spec(:,14)/(max(max(Osram_Led_Spec(:,14))));
Spec_PCLime = Osram_Led_Spec(:,15)/(max(max(Osram_Led_Spec(:,15))));
Spec_Green_517nm = Osram_Led_Spec(:,16)/(max(max(Osram_Led_Spec(:,16))));
Spec_Amber_593nm = Osram_Led_Spec(:,17)/(max(max(Osram_Led_Spec(:,17))));
Spec_PCAmber = Osram_Led_Spec(:,18)/(max(max(Osram_Led_Spec(:,18))));
Spec_Red_633nm = Osram_Led_Spec(:,19)/(max(max(Osram_Led_Spec(:,19))));
Spec_DeepRed_660nm = Osram_Led_Spec(:,20)/(max(max(Osram_Led_Spec(:,20))));
Spec_FarRed_740nm = Osram_Led_Spec(:,21)/(max(max(Osram_Led_Spec(:,21))));
Spec_InfraRed_850nm = Osram_Led_Spec(:,22)/(max(max(Osram_Led_Spec(:,22))));
Spec_DeepInfraRed_940nm = Osram_Led_Spec(:,23)/(max(max(Osram_Led_Spec(:,23))));


figure(1)
hold on
% plot(Osram_lambda,Spec_UVA_365nm,'LineWidth',2, 'Color', '#DAA8FF')
plot(Osram_lambda,Spec_Violet_385nm,'LineWidth',2,'Color','#BE6CFB')
plot(Osram_lambda,Spec_Violet_395nm,'LineWidth',2,'Color','#A62EFF')
plot(Osram_lambda,Spec_Violet_405nm,'LineWidth',2,'Color','#9D00DC')
plot(Osram_lambda,Spec_DeepBlue_436nm,'LineWidth',2,'Color','#001787')
plot(Osram_lambda,Spec_RoyalBlue_453nm,'LineWidth',2,'Color','#0022C6')
plot(Osram_lambda,Spec_DentalBlue_460nm,'LineWidth',2,'Color','#657FFF')
plot(Osram_lambda,Spec_Cyan_500nm,'LineWidth',2,'Color','#70FFFF')
% plot(Osram_lambda,Spec_PCLime,'LineWidth',2,'Color','#68FFC6')
plot(Osram_lambda,Spec_Green_517nm,'LineWidth',2,'Color','#0AE700')
plot(Osram_lambda,Spec_Amber_593nm,'LineWidth',2,'Color','#FFBF00')
% plot(Osram_lambda,Spec_PCAmber,'LineWidth',2,'Color','#FF9E43')
plot(Osram_lambda,Spec_Red_633nm,'LineWidth',2,'Color','#FF2B2B')
plot(Osram_lambda,Spec_DeepRed_660nm,'LineWidth',2,'Color','#E90000')
plot(Osram_lambda,Spec_FarRed_740nm,'LineWidth',2,'Color','#A60000')
plot(Osram_lambda,Spec_InfraRed_850nm,'LineWidth',2,'Color','#5B1B1B')
xlim([0.35 1])
ylim([0 1])
xlabel('Wavelength (\mum)')
title('OSRAM LEDs (Colored)')
legend('Violet (385 nm)','Violet (395 nm)','Violet (405 nm)','Deep Blue (436 nm)', ...
    'Royal Blue (453 nm)','Dental Blue (460 nm)','Cyan (500 nm)','Green (517 nm)','Amber (593 nm)', ...
    'Red (633 nm)','Deep Red (660 nm)','Far Red (740 nm)','Infrared (850 nm)')
%%
figure(2)
hold on
plot(Osram_lambda,Spec_Warm_White_2200K,'LineWidth',2, 'Color' , '#FFC372' )
plot(Osram_lambda,Spec_Warm_White_3000K,'LineWidth',2, 'Color' , '#FFD7A2')
plot(Osram_lambda,Spec_Neutral_White_4000K,'LineWidth',2, 'Color' , '#ECEAE7')
plot(Osram_lambda,Spec_Cool_White_5500K,'LineWidth',2, 'Color', '#D6E8E8')
plot(Osram_lambda,Spec_Cool_White_6500K,'LineWidth',2, 'Color', '#DEF8F8')
xlim([0.35 0.8])
ylim([0 1])
xlabel('Wavelength (\mum)')
title('OSRAM LEDs (White)')
legend('Warm White 2200K','Warm White 3000K','Neutral White 4000K','Cool White 5500K','Cool White 6500K')
%%

%Reflectivity spectrum for colored LEDs at 0 nm

for i = 1:2691
Ref_Spec_WW_2200K = Gamma(i,valtoindex_L(0))'.*Spec_Warm_White_2200K; %Reflectivity spectrum for red at 0 nm (R_r) 
Ref_Spec_WW_3000K = Gamma(i,valtoindex_L(0))'.*Spec_Warm_White_3000K; 
Ref_Spec_NW_4000K = Gamma(i,valtoindex_L(0))'.*Spec_Neutral_White_4000K;  
Ref_Spec_CW_5500K = Gamma(i,valtoindex_L(0))'.*Spec_Cool_White_5500K; 
Ref_Spec_CW_6500K = Gamma(i,valtoindex_L(0))'.*Spec_Cool_White_6500K;  
Ref_Spec_UVA_365nm = Gamma(i,valtoindex_L(0))'.*Spec_UVA_365nm;  
Ref_Spec_Violet_385nm = Gamma(i,valtoindex_L(0))'.*Spec_Violet_385nm;  
Ref_Spec_Violet_395nm = Gamma(i,valtoindex_L(0))'.*Spec_Violet_395nm; 
Ref_Spec_Violet_405nm = Gamma(i,valtoindex_L(0))'.*Spec_Violet_405nm;  
Ref_Spec_DeepBlue_436nm = Gamma(i,valtoindex_L(0))'.*Spec_DeepBlue_436nm;  
Ref_Spec_RoyalBlue_453nm = Gamma(i,valtoindex_L(0))'.*Spec_RoyalBlue_453nm;  
Ref_Spec_DentalBlue_460nm = Gamma(i,valtoindex_L(0))'.*Spec_DentalBlue_460nm;  
Ref_Spec_Cyan_500nm = Gamma(i,valtoindex_L(0))'.*Spec_Cyan_500nm;  
Ref_Spec_PCLime = Gamma(i,valtoindex_L(0))'.*Spec_PCLime; 
Ref_Spec_Green_517nm = Gamma(i,valtoindex_L(0))'.*Spec_Green_517nm; 
Ref_Spec_Amber_593nm = Gamma(i,valtoindex_L(0))'.*Spec_Amber_593nm; 
Ref_Spec_PCAmber = Gamma(i,valtoindex_L(0))'.*Spec_PCAmber;
Ref_Spec_Red_633nm = Gamma(i,valtoindex_L(0))'.*Spec_Red_633nm; 
Ref_Spec_DeepRed_660nm = Gamma(i,valtoindex_L(0))'.*Spec_DeepRed_660nm;
Ref_Spec_FarRed_740nm = Gamma(i,valtoindex_L(0))'.*Spec_FarRed_740nm; 
Ref_Spec_InfraRed_850nm = Gamma(i,valtoindex_L(0))'.*Spec_InfraRed_850nm; 
Ref_Spec_DeepInfraRed_940nm = Gamma(i,valtoindex_L(0))'.*Spec_DeepInfraRed_940nm; 
end 

Ref_Spec_WW_2200K = Ref_Spec_WW_2200K(1:2691);
Ref_Spec_WW_3000K = Ref_Spec_WW_3000K(1:2691);
Ref_Spec_NW_4000K = Ref_Spec_NW_4000K(1:2691);
Ref_Spec_CW_5500K = Ref_Spec_CW_5500K(1:2691);
Ref_Spec_CW_6500K = Ref_Spec_CW_6500K(1:2691);
Ref_Spec_UVA_365nm = Ref_Spec_UVA_365nm(1:2691);
Ref_Spec_Violet_385nm = Ref_Spec_Violet_385nm(1:2691);
Ref_Spec_Violet_395nm = Ref_Spec_Violet_395nm(1:2691);
Ref_Spec_Violet_405nm = Ref_Spec_Violet_405nm(1:2691);
Ref_Spec_DeepBlue_436nm = Ref_Spec_DeepBlue_436nm(1:2691);
Ref_Spec_RoyalBlue_453nm = Ref_Spec_RoyalBlue_453nm(1:2691);
Ref_Spec_DentalBlue_460nm = Ref_Spec_DentalBlue_460nm(1:2691);
Ref_Spec_Cyan_500nm = Ref_Spec_Cyan_500nm(1:2691);
Ref_Spec_PCLime = Ref_Spec_PCLime(1:2691);
Ref_Spec_Green_517nm = Ref_Spec_Green_517nm(1:2691);
Ref_Spec_Amber_593nm = Ref_Spec_Amber_593nm(1:2691);
Ref_Spec_PCAmber = Ref_Spec_PCAmber(1:2691);
Ref_Spec_Red_633nm = Ref_Spec_Red_633nm(1:2691);
Ref_Spec_DeepRed_660nm = Ref_Spec_DeepRed_660nm(1:2691);
Ref_Spec_FarRed_740nm = Ref_Spec_FarRed_740nm(1:2691);
Ref_Spec_InfraRed_850nm = Ref_Spec_InfraRed_850nm(1:2691);
Ref_Spec_DeepInfraRed_940nm = Ref_Spec_DeepInfraRed_940nm(1:2691);

%%

I_ref_WW_2200K = [];
I_ref_WW_3000K = []; 
I_ref_NW_4000K = []; 
I_ref_CW_5500K = [];
I_ref_CW_6500K = [];
I_ref_UVA_365nm = [];
I_ref_Violet_385nm = [];
I_ref_Violet_395nm = [];
I_ref_Violet_405nm = [];
I_ref_DeepBlue_436nm = [];
I_ref_RoyalBlue_453nm = [];
I_ref_DentalBlue_460nm = [];
I_ref_Cyan_500nm = [];
I_ref_PCLime = [];
I_ref_Green_517nm = [];
I_ref_Amber_593nm = [];
I_ref_PCAmber = [];
I_ref_Red_633nm = [];
I_ref_DeepRed_660nm = [];
I_ref_FarRed_740nm = [];
I_ref_InfraRed_850nm = [];
I_ref_DeepInfraRed_940nm = [];

%%

 %Multiply reflectance with I_Out for every thickness and sum to find total reflected intensity

for i=1:2691
    
    I_ref_WW_2200K(:,i) = sum(Gamma(:,i).*Ref_Spec_WW_2200K)./100;         
    I_ref_WW_3000K(:,i) = sum(Gamma(:,i).*Ref_Spec_WW_3000K)./100;     
    I_ref_NW_4000K(:,i) = sum(Gamma(:,i).*Ref_Spec_NW_4000K)./100;        
    I_ref_CW_5500K(:,i) = sum(Gamma(:,i).*Ref_Spec_CW_5500K)./100;
    I_ref_CW_6500K(:,i) = sum(Gamma(:,i).*Ref_Spec_CW_6500K)./100;
    I_ref_UVA_365nm(:,i) = sum(Gamma(:,i).*Ref_Spec_UVA_365nm)./100;
    I_ref_Violet_385nm(:,i) = sum(Gamma(:,i).*Ref_Spec_Violet_385nm)./100;
    I_ref_Violet_395nm(:,i) = sum(Gamma(:,i).*Ref_Spec_Violet_395nm)./100;
    I_ref_Violet_405nm(:,i) = sum(Gamma(:,i).*Ref_Spec_Violet_405nm)./100;
    I_ref_DeepBlue_436nm(:,i) = sum(Gamma(:,i).*Ref_Spec_DeepBlue_436nm)./100;
    I_ref_RoyalBlue_453nm(:,i) = sum(Gamma(:,i).*Ref_Spec_RoyalBlue_453nm)./100;
    I_ref_DentalBlue_460nm(:,i) = sum(Gamma(:,i).*Ref_Spec_DentalBlue_460nm)./100;
    I_ref_Cyan_500nm(:,i) = sum(Gamma(:,i).*Ref_Spec_Cyan_500nm)./100;
    I_ref_PCLime(:,i) = sum(Gamma(:,i).*Ref_Spec_PCLime)./100;
    I_ref_Green_517nm(:,i) = sum(Gamma(:,i).*Ref_Spec_Green_517nm)./100;
    I_ref_Amber_593nm(:,i) = sum(Gamma(:,i).*Ref_Spec_Amber_593nm)./100;
    I_ref_PCAmber(:,i) = sum(Gamma(:,i).*Ref_Spec_PCAmber)./100;
    I_ref_Red_633nm(:,i) = sum(Gamma(:,i).*Ref_Spec_Red_633nm)./100;
    I_ref_DeepRed_660nm(:,i) = sum(Gamma(:,i).*Ref_Spec_DeepRed_660nm)./100;
    I_ref_FarRed_740nm(:,i) = sum(Gamma(:,i).*Ref_Spec_FarRed_740nm)./100;
    I_ref_InfraRed_850nm(:,i) = sum(Gamma(:,i).*Ref_Spec_InfraRed_850nm./100);
    I_ref_DeepInfraRed_940nm(:,i) = sum(Gamma(:,i).*Ref_Spec_DeepInfraRed_940nm./100);
    %% 
end 

%%


figure(3)
hold on
% plot(L,I_ref_UVA_365nm,'LineWidth',2, 'Color', '#DAA8FF')
plot(L,I_ref_Violet_385nm,'LineWidth',2,'Color','#BE6CFB')
plot(L,I_ref_Violet_395nm,'LineWidth',2,'Color','#A62EFF')
plot(L,I_ref_Violet_405nm,'LineWidth',2,'Color','#9D00DC')
plot(L,I_ref_DeepBlue_436nm,'LineWidth',2,'Color','#001787')
plot(L,I_ref_RoyalBlue_453nm,'LineWidth',2,'Color','#0022C6')
plot(L,I_ref_DentalBlue_460nm,'LineWidth',2,'Color','#657FFF')
plot(L,I_ref_Cyan_500nm,'LineWidth',2,'Color','#70FFFF')
% plot(L,I_ref_PCLime,'LineWidth',2,'Color','#68FFC6')
plot(L,I_ref_Green_517nm,'LineWidth',2,'Color','#0AE700')
plot(L,I_ref_Amber_593nm,'LineWidth',2,'Color','#FFBF00')
% plot(L,I_ref_PCAmber,'LineWidth',2,'Color','#FF9E43')
plot(L,I_ref_Red_633nm,'LineWidth',2,'Color','#FF2B2B')
plot(L,I_ref_DeepRed_660nm,'LineWidth',2,'Color','#E90000')
plot(L,I_ref_FarRed_740nm,'LineWidth',2,'Color','#A60000')
plot(L,I_ref_InfraRed_850nm,'LineWidth',2,'Color','#5B1B1B')
plot(L,I_ref_DeepInfraRed_940nm,'LineWidth',2,'Color','#5B1B1B')
title('Reflected Intensity (I_O_u_t) for Colored LEDs from 0 nm to 300 nm')
xlabel('L (\mum)','FontSize',16);
ylim([0 0.16])
ylabel('Reflected Intensity (I_O_u_t)','FontSize',16);
legend('Violet (385 nm)','Violet (395 nm)','Violet (405 nm)','Deep Blue (436 nm)', ...
    'Royal Blue (453 nm)','Dental Blue (460 nm)','Cyan (500 nm)','Green (517 nm)','Amber (593 nm)', ...
    'Red (633 nm)','Deep Red (660 nm)','Far Red (740 nm)','Infrared (850 nm)')
hold off

figure(4)
hold on
plot(L,I_ref_WW_2200K,'LineWidth',2, 'Color' , '#FFC372')
plot(L,I_ref_WW_3000K,'LineWidth',2, 'Color' , '#FFD7A2')
plot(L,I_ref_NW_4000K,'LineWidth',2, 'Color' , '#ECEAE7')   
plot(L,I_ref_CW_5500K,'LineWidth',2, 'Color', '#D6E8E8')
plot(L,I_ref_CW_6500K,'LineWidth',2, 'Color', '#DEF8F8')
title('Reflected Intensity (I_O_u_t) for White LEDs from 0 nm to 300 nm')
ylim([0 0.6])
xlabel('L (\mum)','FontSize',16);
ylabel('Reflected Intensity (I_O_u_t)','FontSize',16);
legend('Warm White 2200K','Warm White 3000K','Neutral White 4000K','Cool White 5500K','Cool White 6500K')
hold off


% figure(5)
% hold on
% plot(L,I_ref_Red_633nm,'r','LineWidth',2)
% plot(L,I_ref_DentalBlue_460nm,'b','LineWidth',2)
% plot(L,I_ref_Green_517nm,'g','LineWidth',2)
% xline(0.12)
% title('Red, Green, and Blue (Current LEDs)')
% legend('Red (633 nm)','Dental Blue (460 nm)','Green (517 nm)','location','northeast')

% figure(6)
% hold on
% plot(L,I_ref_Violet_405nm,'LineWidth',2,'Color','#9D00DC')
% plot(L,I_ref_FarRed_740nm,'LineWidth',2,'Color','#A60000')
% plot(L,I_ref_InfraRed_850nm,'LineWidth',2,'Color','#5B1B1B')
% xline(0.12)
% title('Violet, Far Red, and Infra Red (Good at 120 nm)')
% legend('Violet (405 nm)','Far Red (740 nm)','Infra Red (850 nm)','location','northeast')
% 
figure(7)
hold on
plot(L,I_ref_DentalBlue_460nm,'b','LineWidth',2)
plot(L,I_ref_Red_633nm,'LineWidth',2,'Color','#FF2B2B')
plot(L,I_ref_DeepInfraRed_940nm,'LineWidth',2,'Color','#5B1B1B')
xline(0.12)
title('dental blue, Far Red, and Infra Red (Good at 120 nm)')
legend('Dental Blue (460 nm)','Far Red (740 nm)','Infra Red (850 nm)','location','northeast')

load("Bulk_Effect_Data.mat")

toc


