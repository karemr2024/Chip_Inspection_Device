clc; clear all; close all;
load('Osram_Data.mat')

Gamma1_n_1 = [];
Z1_n_1 = [];
Gamma1_n_1_33 = [];
Z1_n_1_33 = [];
Gamma1_n_1_34 = [];
Z1_n_1_34 = [];

for i = 1:numel(lambda)
for j = 1:numel(L)
[Gamma1_n_1(i,j),Z1_n_1(i,j)] = multidiels([1;nrSiO(1,i);nrSi(1,i)],L(j).*nrSiO(1,i),lambda(1,i));    
[Gamma1_n_1_33(i,j),Z1_n_1_33(i,j)] = multidiels([1.33;nrSiO(1,i);nrSi(1,i)],L(j).*nrSiO(1,i),lambda(1,i));
[Gamma1_n_1_34(i,j),Z1_n_1_34(i,j)] = multidiels([1.34;nrSiO(1,i);nrSi(1,i)],L(j).*nrSiO(1,i),lambda(1,i));
end
end
%%
Gamma_n_1 = conj(Gamma1_n_1).*Gamma1_n_1;
Gamma_n_1_33 = conj(Gamma1_n_1_33).*Gamma1_n_1_33; %Multiply Gamma with conjugate to get rid of imaginary component
Gamma_n_1_34 = conj(Gamma1_n_1_34).*Gamma1_n_1_34;

%%

figure(1)
hold on
plot(L,Gamma_n_1(valtoindex_lambda(cw_b),:),'b','LineWidth',2)
plot(L,Gamma_n_1(valtoindex_lambda(cw_g),:),'g','LineWidth',2)
plot(L,Gamma_n_1(valtoindex_lambda(cw_r),:),'r','LineWidth',2)
plot(L,Gamma_n_1_33(valtoindex_lambda(cw_b),:),'b','LineWidth',2)
plot(L,Gamma_n_1_34(valtoindex_lambda(cw_b),:),'b','LineWidth',2)
plot(L,Gamma_n_1_33(valtoindex_lambda(cw_g),:),'g','LineWidth',2)
plot(L,Gamma_n_1_34(valtoindex_lambda(cw_g),:),'g','LineWidth',2)
plot(L,Gamma_n_1_33(valtoindex_lambda(cw_r),:),'r','LineWidth',2)
plot(L,Gamma_n_1_34(valtoindex_lambda(cw_r),:),'r','LineWidth',2)
xlabel("Oxide Thickness (um)")
ylabel('Reflectance')
legend('Blue','Green','Red','Blue, ref index = 1.33')
title("Reflectance at refractive index n = 1")

figure(2)
hold on
plot(L,Gamma_n_1(valtoindex_lambda(cw_b),:),'r','LineWidth',2)
plot(L,Gamma_n_1_33(valtoindex_lambda(cw_b),:),'c','LineWidth',2)
plot(L,Gamma_n_1_34(valtoindex_lambda(cw_b),:),'m','LineWidth',2)
xlabel("Oxide Thickness (um)")
ylabel('Reflectance')
legend('Refractive Index n = 1', 'Refractive index n = 1.33','Refractive Index n = 1.34')
title("Reflectance for different refractive index for blue")

figure(3)
hold on
plot(L,Gamma_n_1(valtoindex_lambda(cw_g),:),'r','LineWidth',2)
plot(L,Gamma_n_1_33(valtoindex_lambda(cw_g),:),'c','LineWidth',2)
plot(L,Gamma_n_1_34(valtoindex_lambda(cw_g),:),'m','LineWidth',2)
xlabel("Oxide Thickness (um)")
ylabel('Reflectance')
legend('Refractive Index n = 1', 'Refractive index n = 1.33','Refractive Index n = 1.34')
title("Reflectance for different refractive index for green")


