function est_L = reftocurve(RefRed, RefGreen, RefBlue)

load("Imaging_Data.mat")
RefRed = [cw_r, RefRed]; %Reflectance Point at Red
RefGreen = [cw_g, RefGreen]; %Reflectance Point at Green
RefBlue = [cw_b, RefBlue]; %Reflectance Point at Blue

%Subtract reflectance value of theoretical value from the measured reflectance at R,G,B 
for i = 1:1370
Diff_at_B(:,i) = abs(RefBlue - Gamma(valtoindex_lambda(cw_b),i)); %karelerini al
Diff_at_G(:,i) = abs(RefGreen - Gamma(valtoindex_lambda(cw_g),i));
Diff_at_R(:,i) = abs(RefRed - Gamma(valtoindex_lambda(cw_r),i));

%Difference Of Measured Reflectance Points to theoretical reflectance
%points at each thickness for R,G,B
L_Diff(i,:) = [Diff_at_B(2,i) Diff_at_G(2,i) Diff_at_R(2,i)]; 
end
sqr_L_Diff = (L_Diff).^2;
lsqr_sum_L_Diff = sum(((L_Diff).^2),2); %average of R,G,B differences at each thickness
est_in = find(lsqr_sum_L_Diff == min(lsqr_sum_L_Diff)); %find index of minimum difference
est_L = L(est_in); %find estimated thickness
end