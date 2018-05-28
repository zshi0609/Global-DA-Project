clear all;
close all;
format long e;
RandSeed=clock;
rand('seed',RandSeed(6));
%% rfactor
rfactor=3; % will test 3 to 10
scalar_hour_year=365*24;
depth_MIMICS=1; % UNIT: meter
scalar_concen_conten=1000;% convert mg C/cm3 to g C/m2 with the depth of 1 meter
Input_litterTOsoil1=load('Input_litterTOsoil.mat');
	sum_litter1TOsoil1_content=Input_litterTOsoil1.sum_litter1TOsoil1_content;
	sum_litter2TOsoil1_content=Input_litterTOsoil1.sum_litter2TOsoil1_content;
	sum_litter3TOsoil2_content=Input_litterTOsoil1.sum_litter3TOsoil2_content;
    total_input_to_soil=sum_litter1TOsoil1_content+sum_litter2TOsoil1_content+sum_litter3TOsoil2_content;
%% identifying some common parameters
% MGE; 
%% input
area_grid1=load('area_grid.mat');
    area_grid=area_grid1.area_grid;
f_clay1=load ('f_clay.mat');
    f_clay=f_clay1.f_clay;
% estimating_f_met (see the subroutine) % this is a function
f_met1=load ('f_met.mat');%  Unit: g C/ m2/yr for litter input
    f_met=f_met1.f_met;
    meanInputToLitter1=f_met1.meanInputToLitter1;
        meanInputToLitter1(meanInputToLitter1<=0)=NaN; % attention: get rid of zero input
    meanInputToLitter2And3=f_met1.meanInputToLitter2And3;
        meanInputToLitter2And3(meanInputToLitter2And3<=0)=NaN;% attention: get rid of zero input
Tsoil1=load ('Tsoil.mat'); % mean annual soil temperature / will be extracted from long-term model output (script: Tsoil)
    Tsoil=Tsoil1.Tsoil;
mask_NCSCD1=load('mask_NCSCD.mat');
    mask_NCSCD=mask_NCSCD1.mask_NCSCD;
mask_mic1=load('mask_mic.mat');
    mask_mic=mask_mic1.mask_mic;
mask_litter_input1=load('mask_litter_input.mat');
    mask_litter_input=mask_litter_input1.mask_litter_input;
%% parameters: start
%parameters explanation
fraction_r=0.5; % fraction of microbial biomass of r-strategy

%Pscalar_a & Pscalar_b;the two coefficients of the function
    Pscalar=(2*exp(-2*sqrt(f_clay))).\1;
    
%f_i_met &f_i_struc: allocation of input to SOMp & SOMc
    f_i_met=0.05; % proportion of litterfall goes into SOMp
    f_i_struc=0.05; % proportion of litterfall goes into SOMc 
    
% Vslope & Vint: slope and intercept for Vmax calculation: Vmax=exp(Vslope*T+Vint)*a_v*Vmod;
% a_V: tuning coefficient for Vmax Vmax=exp(Vslope*T+Vint)*a_v*Vmod;
% Vmod_r & Vmod_K: modify vmax
    Vslope=0.063; % unit: ln(mg Cs/(mg MIC)/h)/degree
    Vint=5.47; % unit: ln(mg Cs/(mg MIC)/h)    
    a_V=8*10^-6;
    Vmod_r_LITm=10;Vmod_r_LITs_SOMc=2;Vmod_r_SOMa=10;
    Vmod_K_LITm=3;Vmod_K_LITs_SOMc=3;Vmod_K_SOMa=2;
    
% Kslope & Kint:slope&intercept for Km calculation Km=exp(Kslope*T+Kint)*a_K*Kmod;
% Wieder et al. 2015 used the same K slope for r&k microbe
% a_K: tuning coefficient for Km
% Kmod_r & Kmod_K: modify Km
    Kslope_LITm=0.017; Kslope_LITs_SOMc=0.027; Kslope_SOMa=0.017; % unit:ln(mg C cm-3)/degree
    Kint=3.19; % unit: ln(mg C/cm3)
    a_K=10;
    Kmod_r_SOMa_coeff=0.25;
    Kmod_K_SOMa_coeff=0.167;
    Kmod_r_LITm=0.125;Kmod_r_LITs_SOMc=0.5;Kmod_r_SOMa=Kmod_r_SOMa_coeff*Pscalar;
    Kmod_K_LITm=0.5;Kmod_K_LITs_SOMc=0.25;Kmod_K_SOMa=Kmod_K_SOMa_coeff*Pscalar;
    

% D: desporption rate from SOMp to SOMa
coeff_a_for_D=1.5*10^(-5)*scalar_hour_year;coeff_b_for_D=-1.5;%coeff_b_for_D does NOT need to multiply this scalar due to its nature of curvature.


%KO: futher modeify Km for oxidation of SOMc
KO_r=4;KO_K=4;
Tau_mod_K=1; 
Tau_mod_r=1; 
Tau_r=2; % turnover time:in year-1
Tau_K=2; % turnover time:in year-1
f_p_r=0.3*exp(1.3*f_clay);
f_c_r=0.1*exp(-3*f_met);
f_p_K=0.2*exp(0.8*f_clay);
f_c_K=0.3*exp(-3*f_met);
f_a_r=1-f_p_r-f_c_r;
f_a_K=1-f_p_K-f_c_K;
%%% parameters: end

%% observations
%%%%%%%%%%%%%%%% observations:soil carbon %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MIC_raw1=load ('MIC.mat'); % load regridded Microbial biomass (original data is from Xu et al., 2013 GEB), % unit: g C/ m2 generated by the script: RegridMicrobialBiomass.mat
    MIC_raw=MIC_raw1.MIC_raw;
% match coordination of soilcarbon_o & MIC
MIC_1=MIC_raw(1:144,:,:);
MIC_2=MIC_raw(145:288,:,:);
MIC=[MIC_2;MIC_1];
% soil carbon 
soilcarbon_o1=load('soilcarbon_o.mat'); % unit: kg C/m2 top and subsoil carbon content (NOT concentration)
	soilcarbon_o_raw=soilcarbon_o1.soilcarbon_o_raw;
soilcarbon_o1=soilcarbon_o_raw(1:144,:,:);
soilcarbon_o2=soilcarbon_o_raw(145:288,:,:);
soilcarbon_o=[soilcarbon_o2;soilcarbon_o1];
soilcarbon_o=soilcarbon_o*1000; % convert to g C/m2
soilcarbon_o(soilcarbon_o>3*10^8)=NaN; % get rid of the very large value 10^36
soilcarbon_o_top=soilcarbon_o(:,:,1);
soilcarbon_o_sub=soilcarbon_o(:,:,2);
soilcarbon_o_total=soilcarbon_o_top+soilcarbon_o_sub;
soilcarbon_o_total(mask_NCSCD==1)=0;
soilcarbon_o_NCSCD1=load('soilcarbon_o_NCSCD.mat');
    soilcarbon_o_NCSCD30=soilcarbon_o_NCSCD1.soil_carbon_NCSCD30;
    soilcarbon_o_NCSCD100=soilcarbon_o_NCSCD1.soil_carbon_NCSCD100;
    soilcarbon_o_NCSCD200=soilcarbon_o_NCSCD1.soil_carbon_NCSCD200;
    soilcarbon_o_NCSCD300=soilcarbon_o_NCSCD1.soil_carbon_NCSCD300;
soilcarbon_o_NCSCD=soilcarbon_o_NCSCD100;
soilcarbon_o_total_NCSCD1=cat(3,soilcarbon_o_total,soilcarbon_o_NCSCD);
soilcarbon_o_total_NCSCD=nansum(soilcarbon_o_total_NCSCD1,3);
soilcarbon_o_total_NCSCD(soilcarbon_o_total_NCSCD==0)=NaN;
soilcarbon_o_total_NCSCD(mask_mic==0)=NaN;
soilcarbon_o_total_NCSCD(mask_litter_input==0)=NaN;
count_soilcarbon_o_total_NCSCD=sum(soilcarbon_o_total_NCSCD(:)>0);
soilcarbon_o_total_NCSCD_area=soilcarbon_o_total_NCSCD.*area_grid*10^6;
soilcarbon_o_total_NCSCD_area_global=nansum(soilcarbon_o_total_NCSCD_area(:));
% soilcarbon_o_total(mask_NCSCD==1)=NaN;
% soilcarbon_o_NCSCD_tot(mask_NCSCD==0)=0;
% soilcarbon_o_total_NCSCD=soilcarbon_o_total+soilcarbon_o_NCSCD;
log_soilcarbon_o_total_NCSCD=log10(soilcarbon_o_total_NCSCD);
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MCMC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nsimu=500000;                  
upgrade=1;
allow=5;
par=[f_i_met;f_i_struc;Vslope;Vint;Vmod_r_LITs_SOMc;Vmod_r_SOMa;Vmod_K_LITs_SOMc;Vmod_K_SOMa;Kslope_LITs_SOMc;Kslope_SOMa;...
    Kint;Kmod_r_LITs_SOMc;Kmod_r_SOMa_coeff;Kmod_K_LITs_SOMc;Kmod_K_SOMa_coeff;coeff_a_for_D;coeff_b_for_D;Tau_r;Tau_K;fraction_r;KO_r;KO_K];
par_old=par;
Min1=par/rfactor;
Max1=par*rfactor;
Min=Min1;
Max=Max1;
Min(17)=Max1(17);
Max(17)=Min1(17);
Min(20)=0;
Max(20)=1;
diff_par=Max-Min;
J_old = 3000000;
for simu=1:nsimu;
    while (true)
        par_new = par_old+(rand(22,1)-0.5).*diff_par/allow;
        %check if paramaters are within their minima and maxima
            if (par_new(1)>Min(1)&&par_new(1)<Max(1)...
             && par_new(2)>Min(2)&&par_new(2)<Max(2)...
             && par_new(3)>Min(3)&&par_new(3)<Max(3)...
             && par_new(4)>Min(4)&&par_new(4)<Max(4)...
             && par_new(5)>Min(5)&&par_new(5)<Max(5)...
             && par_new(6)>Min(6)&&par_new(6)<Max(6)...
             && par_new(7)>Min(7)&&par_new(7)<Max(7)...
             && par_new(8)>Min(8)&&par_new(8)<Max(8)...
             && par_new(9)>Min(9)&&par_new(9)<Max(9)...
             && par_new(10)>Min(10)&&par_new(10)<Max(10)...
             && par_new(11)>Min(11)&&par_new(11)<Max(11)...
             && par_new(12)>Min(12)&&par_new(12)<Max(12)...
             && par_new(13)>Min(13)&&par_new(13)<Max(13)...
             && par_new(14)>Min(14)&&par_new(14)<Max(14)...
             && par_new(15)>Min(15)&&par_new(15)<Max(15)...
             && par_new(16)>Min(16)&&par_new(16)<Max(16)...
             && par_new(17)>Min(17)&&par_new(17)<Max(17)...
             && par_new(18)>Min(18)&&par_new(18)<Max(18)...
             && par_new(19)>Min(19)&&par_new(19)<Max(19)...
             && par_new(20)>Min(20)&&par_new(20)<Max(20)...
             && par_new(21)>Min(21)&&par_new(21)<Max(21)...
             && par_new(22)>Min(22)&&par_new(22)<Max(22))             
               break;  
            end
    end
    
   %%%%%%%%%%%%%%% parameter space %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    f_i_met=par_new(1);f_i_struc=par_new(2);Vslope=par_new(3);Vint=par_new(4);Vmod_r_LITs_SOMc=par_new(5);    
    Vmod_r_SOMa=par_new(6);Vmod_K_LITs_SOMc=par_new(7);Vmod_K_SOMa=par_new(8);Kslope_LITs_SOMc=par_new(9);Kslope_SOMa=par_new(10);
    Kint=par_new(11);Kmod_r_LITs_SOMc=par_new(12);Kmod_r_SOMa_coeff=par_new(13);Kmod_K_LITs_SOMc=par_new(14);Kmod_K_SOMa_coeff=par_new(15);
    coeff_a_for_D=par_new(16);coeff_b_for_D=par_new(17);Tau_r=par_new(18);Tau_K=par_new(19);fraction_r=par_new(20);KO_r=par_new(21);KO_K=par_new(22);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Kmod_r_SOMa=Kmod_r_SOMa_coeff*Pscalar;
    Kmod_K_SOMa=Kmod_K_SOMa_coeff*Pscalar;
%% The MIMICS model
MICr=MIC*fraction_r;
MICk=MIC*(1-fraction_r);
% D: desporption rate from SOMp to SOMa
D=coeff_a_for_D*exp(coeff_b_for_D*f_clay);

% Parameters for calculating SOMc
Vmax_r2=exp(Vslope*Tsoil+Vint)*a_V*Vmod_r_LITs_SOMc;% note & flag: Vmax_r2 is for both uptake of LITs && SOMc (Wieder 2015 GMD page1800,last paragraph) unit: (mg Cs/mg MIC/h);
Km_r2=exp(Kslope_LITs_SOMc*Tsoil+Kint)*a_K*Kmod_r_LITs_SOMc;%unit: (mg C /cm3)
Vmax_K2=exp(Vslope*Tsoil+Vint)*a_V*Vmod_K_LITs_SOMc;%unit: (mg Cs/mg MIC/h);
Km_K2=exp(Kslope_LITs_SOMc*Tsoil+Kint)*a_K*Kmod_K_LITs_SOMc;%unit: (mg C /cm3)


% Parameters for calculating SOMa
Vmax_r3=exp(Vslope*Tsoil+Vint)*a_V*Vmod_r_SOMa;%unit: (mg Cs/mg MIC/h);
Km_r3=(exp(Kslope_SOMa*Tsoil+Kint)*a_K).*Kmod_r_SOMa;%unit: (mg C /cm3)
Vmax_K3=exp(Vslope*Tsoil+Vint)*a_V*Vmod_K_SOMa;%unit: (mg Cs/mg MIC/h);
Km_K3=(exp(Kslope_SOMa*Tsoil+Kint)*a_K).*Kmod_K_SOMa;%unit: (mg C /cm3)
Decay_MICr_SOM=MICr*Tau_r;                       % eq 4
Decay_MICk_SOM=MICk*Tau_K;                       % eq 8


    
%% calculating steady states
SOMp=(total_input_to_soil*f_i_met+f_p_r.*Decay_MICr_SOM+f_p_K.*Decay_MICk_SOM)./D;
Decay_SOMp_SOMa=SOMp.*D;                          % eq 9

% solve second order equation to derive SS_SOMc & SS_SOMa
    % define all the coefficients for SS_SOMc
    InputToSOMc=total_input_to_soil*f_i_struc+f_c_r.*Decay_MICr_SOM+f_c_K.*Decay_MICk_SOM; %gC/m2/yr
    coeff_a=MICr.*Vmax_r2*scalar_hour_year;% gCs/h*365*24 to g Cs/yr
    coeff_b=KO_r*Km_r2*scalar_concen_conten;% mg C/cm3*1000=g C/ m2
    coeff_c=MICk.*Vmax_K2*scalar_hour_year;% gCs/h*365*24 to g Cs/yr
    coeff_d=KO_K*Km_K2*scalar_concen_conten;% mg C/cm3*1000=g C/ m2
    coeff_e=coeff_a+coeff_c-InputToSOMc;% the coefficient for X square
    coeff_f=coeff_a.*coeff_d+coeff_b.*coeff_c-(coeff_b+coeff_d).*InputToSOMc;% the coefficient for X 
    coeff_g=-(coeff_b.*coeff_d).*InputToSOMc;% the constant 
    
  SOMc=(-coeff_f+sqrt(coeff_f.^2-4*coeff_e.*coeff_g))./(2*coeff_e);
  SOMc(imag(SOMc)~=0) = NaN;% Complex numbers
  SOMc(SOMc<0)=NaN;% negative numbers
Uptake_SOMc_SOMa=((MICr.*(Vmax_r2*scalar_hour_year)).*SOMc)./(KO_r*Km_r2*scalar_concen_conten+SOMc)+((MICk.*(Vmax_K2*scalar_hour_year)).*SOMc)./(KO_K*Km_K2*scalar_concen_conten+SOMc); % eq 10 with unit conversion scalars to g C/m2/yr
	% define all the coefficients for SS_SOMa
    InputToSOMa=f_a_r.*Decay_MICr_SOM+f_a_K.*Decay_MICk_SOM+Decay_SOMp_SOMa+Uptake_SOMc_SOMa;
    coeff_h=MICr.*Vmax_r3*scalar_hour_year;
    coeff_i=Km_r3*scalar_concen_conten;
    coeff_j=MICk.*Vmax_K3*scalar_hour_year;
    coeff_k=Km_K3*scalar_concen_conten;
    coeff_l=coeff_h+coeff_j-InputToSOMa;% the coefficient for X square
    coeff_m=coeff_h.*coeff_k+coeff_i.*coeff_j-(coeff_i+coeff_k).*InputToSOMa;% the coefficient for X 
    coeff_n=(-coeff_i.*coeff_k).*InputToSOMa;% the constant 
  
    SOMa=(-coeff_m+sqrt(coeff_m.^2-4*coeff_l.*coeff_n))./(2*coeff_l);
  SOMa(imag(SOMa)~=0) = NaN;% Complex numbers
  SOMa(SOMa<0)=NaN;% negative numbers

Soilcarbon_M=SOMp+SOMc+SOMa;
Soilcarbon_M(Soilcarbon_M<=0)=NaN;
Soilcarbon_M(mask_mic==0)=NaN;
Soilcarbon_M(mask_litter_input==0)=NaN;
Soilcarbon_M_area=Soilcarbon_M.*area_grid*10^6;
SS=nansum(Soilcarbon_M_area(:));

logged_soilcarbon_M=log10(Soilcarbon_M);
scalar_sd=0.5;
diff=logged_soilcarbon_M-log_soilcarbon_o_total_NCSCD;
J1=(diff.^2)./(2*(scalar_sd*log_soilcarbon_o_total_NCSCD).^2);
J_new=nansum(J1(:));
delta_J = J_new-J_old;
check_point1=sum(SOMa(:)>0);   
check_point2=sum(SOMc(:)>0);  
if min(1,exp(-delta_J))>rand && sum(SOMa(:)>0)>9000
   Parameters_keep(:,upgrade)=par_new;
   J_keep(upgrade)=J_new;
   count_SOMc(upgrade)=sum(SOMc(:)>0);
   count_SOMa(upgrade)=sum(SOMa(:)>0);
   count_SOMp(upgrade)=sum(SOMp(:)>0);
   upgrade=upgrade+1
   par_old=par_new;
   J_old=J_new;  
   ss_record(upgrade)=SS;  
end
 
