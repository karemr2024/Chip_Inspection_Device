clc; clearvars; close all;
load("reftocurve_data.mat")                      %Load Gamma Matrix (Gamma), Thickness (L), Wavelength
                                                 %(lambda), center wavelengths (cw_b, cw_g, cw_r)

Ref_Blue = input("Enter blue reflectance: ");    %Reflectance of Oxide at Blue
Ref_Green = input("Enter green reflectance: ");  %Reflectance of Oxide at Green
Ref_Red = input("Enter red reflectance: ");      %Reflectance of Oxide at Red




esti_L = reftocurve(Ref_Red,Ref_Green,Ref_Blue); %Estimated thickness calculated using 
                                                 %reftocurve function and R,G,B reflectances
                                                 %as inputs

figure(1)

hold on
plot(lambda,Gamma(:,valtoindex_L(abs(esti_L))),'m--','LineWidth',2)           %Estimated nm thickness
plot(lambda,Gamma(:,valtoindex_L(abs(esti_L+0.001))),'c--','LineWidth',2)     %Estimated +1 nm thickness
plot(lambda,Gamma(:,valtoindex_L(abs(esti_L-0.001))),'y--','LineWidth',2)     %Estimated -1 nm thickness

plot(cw_r,Ref_Red,'r.','MarkerSize',30)
plot(cw_g,Ref_Green,'g.','MarkerSize',30)
plot(cw_b,Ref_Blue,'b.','MarkerSize',30)
xlabel('lambda (\mum)')
ylabel('Reflectance')
xlim([0.4 0.68])
ylim([0 0.45])
legend1 = sprintf('L = %0.2f nm', 1000*esti_L);
legend2 = sprintf('L = %0.2f nm', 1000*(esti_L+0.001));
legend3 = sprintf('L = %0.2f nm', 1000*(esti_L-0.001));
legend(legend1, legend2, legend3,'Ref at R','Ref at G','Ref at B','location','bestoutside')
fprintf('Estimated thickness is %f nm \n',esti_L*1000)