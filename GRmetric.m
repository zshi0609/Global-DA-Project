clear
% clm 4
clm4_chain1=csvread('d:\oco2_NC\clm4_NC_final\accepted_parameters\clm4_chain1.csv');
clm4_chain2=csvread('d:\oco2_NC\clm4_NC_final\accepted_parameters\clm4_chain2.csv');
clm4_chain3=csvread('d:\oco2_NC\clm4_NC_final\accepted_parameters\clm4_chain3.csv');
clm4_chain4=csvread('d:\oco2_NC\clm4_NC_final\accepted_parameters\clm4_chain4.csv');
clm4_chain5=csvread('d:\oco2_NC\clm4_NC_final\accepted_parameters\clm4_chain5.csv');
% take the last 10,000 parameters from each of the chain
aaa=9999;
% parameters_clm4=[clm4_chain1(length(clm4_chain1(:,1))-aaa:length(clm4_chain1(:,1)),:);clm4_chain2(length(clm4_chain2(:,1))-aaa:length(clm4_chain2(:,1)),:);...
%     clm4_chain3(length(clm4_chain3(:,1))-aaa:length(clm4_chain3(:,1)),:);clm4_chain4(length(clm4_chain4(:,1))-aaa:length(clm4_chain4(:,1)),:);...
%     clm4_chain5(length(clm4_chain5(:,1))-aaa:length(clm4_chain5(:,1)),:)];
x1=clm4_chain1(length(clm4_chain1(:,1))-aaa:length(clm4_chain1(:,1)),:);
x2=clm4_chain2(length(clm4_chain2(:,1))-aaa:length(clm4_chain2(:,1)),:);
x3=clm4_chain3(length(clm4_chain3(:,1))-aaa:length(clm4_chain3(:,1)),:);
x4=clm4_chain4(length(clm4_chain4(:,1))-aaa:length(clm4_chain4(:,1)),:);
x5=clm4_chain5(length(clm4_chain5(:,1))-aaa:length(clm4_chain5(:,1)),:);
K=5;N=10000;
x1_bar=mean(x1);x2_bar=mean(x2);x3_bar=mean(x3);x4_bar=mean(x4);x5_bar=mean(x5);
x_bar=(x1_bar+x2_bar+x3_bar+x4_bar+x5_bar)/K;
for par_no=1:10;
sumsquare_x1=(x1(:,par_no)-x1_bar(par_no))'*(x1(:,par_no)-x1_bar(par_no));
sumsquare_x2=(x2(:,par_no)-x2_bar(par_no))'*(x2(:,par_no)-x2_bar(par_no));
sumsquare_x3=(x3(:,par_no)-x3_bar(par_no))'*(x3(:,par_no)-x3_bar(par_no));
sumsquare_x4=(x4(:,par_no)-x4_bar(par_no))'*(x4(:,par_no)-x4_bar(par_no));
sumsquare_x5=(x5(:,par_no)-x5_bar(par_no))'*(x5(:,par_no)-x5_bar(par_no));
Bi=N/(K-1)*((x1_bar(par_no)-x_bar(par_no))^2+(x2_bar(par_no)-x_bar(par_no))^2+(x3_bar(par_no)-x_bar(par_no))^2+(x4_bar(par_no)-x_bar(par_no))^2+(x5_bar(par_no)-x_bar(par_no))^2);
Wi=1/(K*(N-1))*(sumsquare_x1+sumsquare_x2+sumsquare_x3+sumsquare_x4+sumsquare_x5);
GR=sqrt((Wi*(N-1)/N+Bi/N)/Wi);
GR_record(par_no)=GR;
end
%
clear
% clm 45
clm45_chain1=csvread('d:\oco2_NC\clm45_NC_final\accepted_parameters\mcmc_clm45_chain1.csv');
clm45_chain2=csvread('d:\oco2_NC\clm45_NC_final\accepted_parameters\mcmc_clm45_chain2.csv');
clm45_chain3=csvread('d:\oco2_NC\clm45_NC_final\accepted_parameters\mcmc_clm45_chain3.csv');
clm45_chain4=csvread('d:\oco2_NC\clm45_NC_final\accepted_parameters\mcmc_clm45_chain4.csv');
clm45_chain5=csvread('d:\oco2_NC\clm45_NC_final\accepted_parameters\mcmc_clm45_chain5.csv');
% take the last 10,000 parameters from each of the chain
aaa=9999;
% parameters_clm4=[clm4_chain1(length(clm4_chain1(:,1))-aaa:length(clm4_chain1(:,1)),:);clm4_chain2(length(clm4_chain2(:,1))-aaa:length(clm4_chain2(:,1)),:);...
%     clm4_chain3(length(clm4_chain3(:,1))-aaa:length(clm4_chain3(:,1)),:);clm4_chain4(length(clm4_chain4(:,1))-aaa:length(clm4_chain4(:,1)),:);...
%     clm4_chain5(length(clm4_chain5(:,1))-aaa:length(clm4_chain5(:,1)),:)];
x1=clm45_chain1(length(clm45_chain1(:,1))-aaa:length(clm45_chain1(:,1)),:);
x2=clm45_chain2(length(clm45_chain2(:,1))-aaa:length(clm45_chain2(:,1)),:);
x3=clm45_chain3(length(clm45_chain3(:,1))-aaa:length(clm45_chain3(:,1)),:);
x4=clm45_chain4(length(clm45_chain4(:,1))-aaa:length(clm45_chain4(:,1)),:);
x5=clm45_chain5(length(clm45_chain5(:,1))-aaa:length(clm45_chain5(:,1)),:);
K=5;N=10000;
x1_bar=mean(x1);x2_bar=mean(x2);x3_bar=mean(x3);x4_bar=mean(x4);x5_bar=mean(x5);
x_bar=(x1_bar+x2_bar+x3_bar+x4_bar+x5_bar)/K;
for par_no=1:13;
sumsquare_x1=(x1(:,par_no)-x1_bar(par_no))'*(x1(:,par_no)-x1_bar(par_no));
sumsquare_x2=(x2(:,par_no)-x2_bar(par_no))'*(x2(:,par_no)-x2_bar(par_no));
sumsquare_x3=(x3(:,par_no)-x3_bar(par_no))'*(x3(:,par_no)-x3_bar(par_no));
sumsquare_x4=(x4(:,par_no)-x4_bar(par_no))'*(x4(:,par_no)-x4_bar(par_no));
sumsquare_x5=(x5(:,par_no)-x5_bar(par_no))'*(x5(:,par_no)-x5_bar(par_no));
Bi=N/(K-1)*((x1_bar(par_no)-x_bar(par_no))^2+(x2_bar(par_no)-x_bar(par_no))^2+(x3_bar(par_no)-x_bar(par_no))^2+(x4_bar(par_no)-x_bar(par_no))^2+(x5_bar(par_no)-x_bar(par_no))^2);
Wi=1/(K*(N-1))*(sumsquare_x1+sumsquare_x2+sumsquare_x3+sumsquare_x4+sumsquare_x5);
GR=sqrt((Wi*(N-1)/N+Bi/N)/Wi);
GR_record(par_no)=GR;
end
%
clear
% MIMICS
mimics_chain1=csvread('d:\oco2_NC\mimics_NC_final\accepted_parameters\mcmc_mimics_chain1.csv');
mimics_chain2=csvread('d:\oco2_NC\mimics_NC_final\accepted_parameters\mcmc_mimics_chain2.csv');
mimics_chain3=csvread('d:\oco2_NC\mimics_NC_final\accepted_parameters\mcmc_mimics_chain3.csv');
mimics_chain4=csvread('d:\oco2_NC\mimics_NC_final\accepted_parameters\mcmc_mimics_chain4.csv');
mimics_chain5=csvread('d:\oco2_NC\mimics_NC_final\accepted_parameters\mcmc_mimics_chain5.csv');
% take the last 10,000 parameters from each of the chain
aaa=9999;
% parameters_clm4=[clm4_chain1(length(clm4_chain1(:,1))-aaa:length(clm4_chain1(:,1)),:);clm4_chain2(length(clm4_chain2(:,1))-aaa:length(clm4_chain2(:,1)),:);...
%     clm4_chain3(length(clm4_chain3(:,1))-aaa:length(clm4_chain3(:,1)),:);clm4_chain4(length(clm4_chain4(:,1))-aaa:length(clm4_chain4(:,1)),:);...
%     clm4_chain5(length(clm4_chain5(:,1))-aaa:length(clm4_chain5(:,1)),:)];
x1=mimics_chain1(length(mimics_chain1(:,1))-aaa:length(mimics_chain1(:,1)),:);
x2=mimics_chain2(length(mimics_chain2(:,1))-aaa:length(mimics_chain2(:,1)),:);
x3=mimics_chain3(length(mimics_chain3(:,1))-aaa:length(mimics_chain3(:,1)),:);
x4=mimics_chain4(length(mimics_chain4(:,1))-aaa:length(mimics_chain4(:,1)),:);
x5=mimics_chain5(length(mimics_chain5(:,1))-aaa:length(mimics_chain5(:,1)),:);
K=5;N=10000;
x1_bar=mean(x1);x2_bar=mean(x2);x3_bar=mean(x3);x4_bar=mean(x4);x5_bar=mean(x5);
x_bar=(x1_bar+x2_bar+x3_bar+x4_bar+x5_bar)/K;
for par_no=1:22;
sumsquare_x1=(x1(:,par_no)-x1_bar(par_no))'*(x1(:,par_no)-x1_bar(par_no));
sumsquare_x2=(x2(:,par_no)-x2_bar(par_no))'*(x2(:,par_no)-x2_bar(par_no));
sumsquare_x3=(x3(:,par_no)-x3_bar(par_no))'*(x3(:,par_no)-x3_bar(par_no));
sumsquare_x4=(x4(:,par_no)-x4_bar(par_no))'*(x4(:,par_no)-x4_bar(par_no));
sumsquare_x5=(x5(:,par_no)-x5_bar(par_no))'*(x5(:,par_no)-x5_bar(par_no));
Bi=N/(K-1)*((x1_bar(par_no)-x_bar(par_no))^2+(x2_bar(par_no)-x_bar(par_no))^2+(x3_bar(par_no)-x_bar(par_no))^2+(x4_bar(par_no)-x_bar(par_no))^2+(x5_bar(par_no)-x_bar(par_no))^2);
Wi=1/(K*(N-1))*(sumsquare_x1+sumsquare_x2+sumsquare_x3+sumsquare_x4+sumsquare_x5);
GR=sqrt((Wi*(N-1)/N+Bi/N)/Wi);
GR_record(par_no)=GR;
end