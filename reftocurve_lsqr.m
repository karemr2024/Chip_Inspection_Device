function est_L = reftocurve_lsqr(RefRed, RefGreen, RefBlue)

load("Simulation_Data.mat")

%% Get reflectance points at RGB

RefRed = [cw_r, RefRed]; %Reflectance Point at Red
RefGreen = [cw_g, RefGreen]; %Reflectance Point at Green
RefBlue = [cw_b, RefBlue]; %Reflectance Point at Blue

%% Subtract reflectance value of theoretical value from the measured reflectance at RGB 
%% Creates a n-by-3 (3 = RGB) matrix with difference of measured and theoritical reflectance values. 

for i = 1:length(Gamma)
    Diff_at_B(:,i) = abs(RefBlue - Gamma_B(valtoindex_lambda(cw_b),i));
    Diff_at_G(:,i) = abs(RefGreen - Gamma_G(valtoindex_lambda(cw_g),i));
    Diff_at_R(:,i) = abs(RefRed - Gamma_R(valtoindex_lambda(cw_r),i));

    % Diff_at_B(:,i) = abs(RefBlue - Gamma(valtoindex_lambda(cw_b),i));
    % Diff_at_G(:,i) = abs(RefGreen - Gamma(valtoindex_lambda(cw_g),i));
    % Diff_at_R(:,i) = abs(RefRed - Gamma(valtoindex_lambda(cw_r),i));

    L_Diff(i,:) = [Diff_at_B(2,i) Diff_at_G(2,i) Diff_at_R(2,i)]; 

end

lsqr_sum_L_Diff = sum(((L_Diff).^2),2); %average of R,G,B differences at each thickness
est_in = find(lsqr_sum_L_Diff == min(lsqr_sum_L_Diff)); %find index of minimum difference
est_L = L(est_in); %find estimated thickness

end