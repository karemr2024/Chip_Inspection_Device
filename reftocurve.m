function est_L = reftocurve(RefRed, RefGreen, RefBlue)

load("reftocurve_data.mat")

RefRed = [cw_r, RefRed]; %Reflectance Point at Red
RefGreen = [cw_g, RefGreen]; %Reflectance Point at Green
RefBlue = [cw_b, RefBlue]; %Reflectance Point at Blue

for i = 1:1370
Diff_at_B(:,i) = abs(RefBlue - Gamma(valtoindex_lambda(cw_b),i));
Diff_at_G(:,i) = abs(RefGreen - Gamma(valtoindex_lambda(cw_g),i));
Diff_at_R(:,i) = abs(RefRed - Gamma(valtoindex_lambda(cw_r),i));
L_Diff(i,:) = [Diff_at_B(2,i) Diff_at_G(2,i) Diff_at_R(2,i)];
end
avg_L_Diff = mean(L_Diff,2);
est_in = find(avg_L_Diff == min(avg_L_Diff));
est_L = L(est_in);
end

