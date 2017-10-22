
%%
clear all;
close all;
format long e;
RandSeed=clock;
rand('seed',RandSeed(6));
Input_litterTOsoil1=load('Input_litterTOsoil.mat');
	sum_litter1TOsoil1_content=Input_litterTOsoil1.sum_litter1TOsoil1_content;
	sum_litter2TOsoil1_content=Input_litterTOsoil1.sum_litter2TOsoil1_content;
	sum_litter3TOsoil2_content=Input_litterTOsoil1.sum_litter3TOsoil2_content;
%% read in scalars
area_grid1=load('area_grid.mat');
    area_grid=area_grid1.area_grid;
scalar1=load('scalar.mat'); 
	mean_scalar_soildecay=scalar1.mean_scalar_soildecay; 
Tsoil1=load ('Tsoil.mat'); 
    Tsoil=Tsoil1.Tsoil;
%% load auxiliary datasets
mean_sand1=load('mean_sand.mat');
	mean_sand=mean_sand1.mean_sand;
mask_NCSCD1=load('mask_NCSCD.mat');
    mask_NCSCD=mask_NCSCD1.mask_NCSCD;
mask_mic1=load('mask_mic.mat');
    mask_mic=mask_mic1.mask_mic;
mask_litter_input1=load('mask_litter_input.mat');
    mask_litter_input=mask_litter_input1.mask_litter_input;
%% %%%%%%%%%%%%%%%%% Parameter ranges %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%define initial parameter values, minima, and maxima; 9 in total
%coefficients t (t1 & t2),transfer coefficients (fs=[f31;
%...;f12;f32;f13], f21=1-t-f31;baseline decay rate:K0=[K0soil1,K0soil2,K0soil3];
% range of ts in the equation t = t1-t2*0.01*(100-CELLSAND_meanovertime);
% t = 0.85-0.68*0.01*(100-CELLSAND_meanovertime);
min_t1=0; max_t1=1;t1=0.85;
min_t2=0; max_t2=1;t2=0.68;
ts=[t1;t2];
min_ts=[min_t1;min_t2];max_ts=[max_t1;max_t2];
% range of transfer coefficients
%fs=[f31;f12;f32;f13]
f31=0.004;min_f31=0;max_f31=0.01;
f12=0.42;min_f12=0.1;max_f12=0.6;
f32=0.03;min_f32=0;max_f32=0.05;
f13=0.45;min_f13=0.3;max_f13=0.7;
fs=[f31;f12;f32;f13];
min_fs=[min_f31;min_f12;min_f32;min_f13];% refer to Shi et al., 2015 Table 2
max_fs=[max_f31;max_f12;max_f32;max_f13];% refer to Shi et al., 2015 Table 2
% range of decay rate K0=[K0soil1,K0soil2,K0soil3];
% refer to Shi et al., 2015 Table 2 & Hararuk et al., 2014 Table 1
factor=1;
min_K0=[1/factor;0.1/factor;0.001/factor];%[1;0.1;0.001];
max_K0=[15*factor;0.5*factor;0.01*factor];
K0=[7.3;0.2;0.0045]; % from clm4.5 codes CNDecompCascadeMod_CENTURY.F90
min_Q10=1;max_Q10=3;Q10=2; % for range see Hararuk et al., 2014 JGR
%% %%%%%%%%%% I %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% litter input to soil
I_l1_s1=sum_litter1TOsoil1_content;
I_l2_s1=sum_litter2TOsoil1_content;
I_l3_s2=sum_litter3TOsoil2_content;
I(:,:,1)=I_l1_s1+I_l2_s1;
I(:,:,2)=I_l3_s2;
I(:,:,3)=0;
%%%%%%%%%%%%%%%% observations:soil carbon %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
soilcarbon_o_NCSCD=soilcarbon_o_NCSCD100+soilcarbon_o_NCSCD200+soilcarbon_o_NCSCD300;
soilcarbon_o_total_NCSCD1=cat(3,soilcarbon_o_total,soilcarbon_o_NCSCD);
soilcarbon_o_total_NCSCD=nansum(soilcarbon_o_total_NCSCD1,3);
soilcarbon_o_total_NCSCD(soilcarbon_o_total_NCSCD==0)=NaN;
soilcarbon_o_total_NCSCD(mask_mic==0)=NaN;
soilcarbon_o_total_NCSCD(mask_litter_input==0)=NaN;
soilcarbon_o_total_NCSCD_area=soilcarbon_o_total_NCSCD.*area_grid*10^6;
soilcarbon_o_total_NCSCD_area_global=nansum(soilcarbon_o_total_NCSCD_area(:));
log_soilcarbon_o_total_NCSCD=log10(soilcarbon_o_total_NCSCD);
%%%%%%%%%%%%%%%%%%%%%%%% step 1      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% MCMC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters_keep= [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]';
nsimu1=1000000;   
upgrade=1;
allow=5;
par=[ts;fs;K0;Q10];
par_old=par;
Min=[min_ts;min_fs;min_K0;min_Q10];
Max=[max_ts;max_fs;max_K0;max_Q10];
diff_par=Max-Min;
J_old = 3000000;
for simu=1:nsimu1
		while (true)
            %propose new sets of parameters
        par_new = par_old+(rand(10,1)-0.5).*diff_par/allow;
        t1=par_new(1);t2=par_new(2);
        t = t1-t2*0.01*(100-mean_sand);
        f31=par_new(3);
        f21=1-t-f31;
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
             && min(t(:))>0 && min(f21(:))>0)             
               break;  
            end              
        end  
%%%%%%%%%%%%%%% parameter space %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fs=par_new(3:6);
K0=par_new(7:9);
Q10=par_new(10);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T_SCALAR=Q10.^((Tsoil-25)/10);
K_s1=(K0(1)*mean_scalar_soildecay).*T_SCALAR;
K_s2=(K0(2)*mean_scalar_soildecay).*T_SCALAR;
K_s3=(K0(3)*mean_scalar_soildecay).*T_SCALAR;
K_s(:,:,1)=K_s1;
K_s(:,:,2)=K_s2;
K_s(:,:,3)=K_s3;
%%%%%%%%%%%%%%% for equilibrium calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
soilcarbon_ss_emulator=nan(288,192,3);
k=1;
for lon=1:288;
    for lat=1:192; 
         if  soilcarbon_o_total_NCSCD(lon,lat)>0 ; % in order to keep input larger than zero; because the mask_overall applied to soilcarbon_o_total
                k=k+1;
                A=[-1 fs(2) fs(4);f21(lon,lat) -1 0;fs(1) fs(3) -1]; % the transfer matrix  
                soilcarbon_ss_emulator(lon,lat,1:3)=(A*diag(squeeze(K_s(lon,lat,:))))\(squeeze(-I(lon,lat,:)));
         end
    end
end

soilcarbon_ss_emulator_total=sum(soilcarbon_ss_emulator(:,:,:),3);
soilcarbon_ss_emulator_total(soilcarbon_ss_emulator_total==0)=NaN;
soilcarbon_ss_emulator_total(mask_mic==0)=NaN;
soilcarbon_ss_emulator_total(mask_litter_input==0)=NaN;
soilcarbon_ss_emulator_total_area=soilcarbon_ss_emulator_total.*area_grid*10^6;
soilcarbon_ss_emulator_total_area_global=nansum(soilcarbon_ss_emulator_total_area(:));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
log_soilcarbon_M_total=log10(soilcarbon_ss_emulator_total);
scalar_sd=0.5;
diff=log_soilcarbon_M_total-log_soilcarbon_o_total_NCSCD;
J=(diff.^2)./(2*(scalar_sd*log_soilcarbon_o_total_NCSCD).^2);
J_new=nansum(nansum(J));
delta_J = J_new-J_old;
   
if min(1,exp(-delta_J))>rand;
   Parameters_keep(:,upgrade)=par_new;
   J_keep(upgrade)=J_new;
   upgrade=upgrade+1;
   par_old=par_new;
   J_old=J_new; 
   count_record(upgrade)=sum(soilcarbon_ss_emulator_total(:)>0);
   ss_record(upgrade)=soilcarbon_ss_emulator_total_area_global;
end
simu
upgrade
Parameters_rec(:,simu)=par_old;      
J_rec(:,simu)=J_old;   
end
