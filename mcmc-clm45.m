
clear all;
close all;
format long e;
RandSeed=clock;
rand('seed',RandSeed(6));
soildepth1=load('soildepth.mat');
	zisoi=soildepth1.zisoi;
	dz=soildepth1.dz;
	zsoi=soildepth1.zsoi;
	dz_node=soildepth1.dz_node;
soildepth_gridded1=load('soildepth_gridded.mat');
    soildepth_gridded=soildepth_gridded1.soildepth_gridded;
Input_litterTOsoil=load('Input_litterTOsoil.mat');
    mean_litter1TOsoil1=Input_litterTOsoil.mean_litter1TOsoil1;
    mean_litter2TOsoil1=Input_litterTOsoil.mean_litter2TOsoil1;
    mean_litter3TOsoil2=Input_litterTOsoil.mean_litter3TOsoil2;
    I=Input_litterTOsoil.I;
mask_perm2_1=load('mask_perm2.mat'); 
    mask_perm2=mask_perm2_1.mask_perm2;
mask_overall1=load('mask_overall.mat');
	mask_overall=mask_overall1.mask_overall;
mask_NCSCD1=load('mask_NCSCD.mat');
    mask_NCSCD=mask_NCSCD1.mask_NCSCD;
%%
scalar1=load('scalar.mat'); 
    scalar_soil_mean=scalar1.Scalar_soildecay;
sand=load('sand.mat');
    CELLSAND_meanovertime=sand.CELLSAND_meanovertime;
diffusivity=load('diffusivity.mat');
    Diffus_meanovertime=diffusivity.Diffus_meanovertime;
Tsoil1=load ('Tsoil.mat');
    Tsoil=Tsoil1.Tsoil;
%%
%%%%%%%%%%%%%%%%%%% Parameter ranges %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%define initial parameter values, minima, and maxima; 12 in total
%diffusivity (D1 (non-permafrost), D2 (permafrost) ),e-folding depth (zt),
%coefficients t (t1 & t2),transfer coefficients (fs=[f31;
%...;f12;f32;f13], f21=1-t-f31;baseline decay rate:K0=[K0soil1,K0soil2,K0soil3];


% range of diffusivity (from Table 3 Koven et al., 2014)
min_D1=0.3*10^-4;% diffusivity minima 
max_D1=16.58*10^-4;%diffusivity maxima
D1=1*10^-4; %unit: m2/year for non-permafrost
min_D2=0.3*10^-4;% diffusivity minima
max_D2=16.58*10^-4;%diffusivity maxima
D2=5*10^-4; %unit: m2/year for permafrost

% rang of e-folding depth
zt=0.5;
min_zt=0;
max_zt=1;

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
factor=1;% can be used to expand the parameter range
min_K0=[1/factor;0.1/factor;0.001/factor];%[1;0.1;0.001];
max_K0=[15*factor;0.5*factor;0.01*factor];
K0=[7.3;0.2;0.0045]; % unit: year from clm4.5 codes CNDecompCascadeMod_CENTURY.F90
min_Q10=1; max_Q10=3;Q10=1.5;



%%%%%%%%%%%%%%%% observations:soil carbon %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% match coordination of soilcarbon_o 
soilcarbon_o1=load('soilcarbon_o.mat'); % unit: kg C/m2 top and subsoil carbon content (NOT concentration)
	soilcarbon_o_raw=soilcarbon_o1.soilcarbon_o_raw;
soilcarbon_o1=soilcarbon_o_raw(1:144,:,:);
soilcarbon_o2=soilcarbon_o_raw(145:288,:,:);
soilcarbon_o=[soilcarbon_o2;soilcarbon_o1];
soilcarbon_o=soilcarbon_o*1000; % convert to g C/m2
soilcarbon_o(soilcarbon_o>3*10^8)=NaN; % get rid of the very large value 10^36
soilcarbon_o_NCSCD1=load('soilcarbon_o_NCSCD.mat');
    soilcarbon_o_NCSCD30=soilcarbon_o_NCSCD1.soil_carbon_NCSCD30;
    soilcarbon_o_NCSCD100=soilcarbon_o_NCSCD1.soil_carbon_NCSCD100;
    soilcarbon_o_NCSCD200=soilcarbon_o_NCSCD1.soil_carbon_NCSCD200;soilcarbon_o_NCSCD200(soilcarbon_o_NCSCD200==0)=NaN;
    soilcarbon_o_NCSCD300=soilcarbon_o_NCSCD1.soil_carbon_NCSCD300;soilcarbon_o_NCSCD300(soilcarbon_o_NCSCD300==0)=NaN;
soilcarbon_o_top=soilcarbon_o(:,:,1);
soilcarbon_o_top(mask_overall(:,:,1)==0)=NaN;
soilcarbon_o_top(mask_NCSCD==1)=0;
soilcarbon_o_top_NCSCD1=cat(3,soilcarbon_o_top,soilcarbon_o_NCSCD30);
soilcarbon_o_top_NCSCD=nansum(soilcarbon_o_top_NCSCD1,3);
soilcarbon_o_top_NCSCD(soilcarbon_o_top_NCSCD==0)=NaN;
soilcarbon_o_total=soilcarbon_o(:,:,1)+soilcarbon_o(:,:,2);
soilcarbon_o_total(mask_NCSCD==1)=0;
soilcarbon_o_total_NCSCD1=cat(3,soilcarbon_o_total,soilcarbon_o_NCSCD100);
soilcarbon_o_total_NCSCD=nansum(soilcarbon_o_total_NCSCD1,3);
soilcarbon_o_total_NCSCD(soilcarbon_o_total_NCSCD==0)=NaN;
soilcarbon_o_top_NCSCD(mask_overall(:,:,1)==0)=NaN;
soilcarbon_o_total_NCSCD(mask_overall(:,:,1)==0)=NaN;
soilcarbon_o_NCSCD200(mask_overall(:,:,1)==0)=NaN;
soilcarbon_o_NCSCD300(mask_overall(:,:,1)==0)=NaN;
log_soilcarbon_o_top_NCSCD=log10(soilcarbon_o_top_NCSCD);
log_soilcarbon_o_total_NCSCD=log10(soilcarbon_o_total_NCSCD);
log_soilcarbon_o_NCSCD200=log10(soilcarbon_o_NCSCD200);
log_soilcarbon_o_NCSCD300=log10(soilcarbon_o_NCSCD300);
%%%%%%%%%%%%%%%%%%%%%%%% step 1      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% MCMC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters_keep= [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]';
nsimu1=1000000;                  
upgrade=1;
allow=5;
par=[D1;D2;zt;ts;fs;K0;Q10];
par_old=par;
Min=[min_D1;min_D2;min_zt;min_ts;min_fs;min_K0;min_Q10];
Max=[max_D1;max_D2;max_zt;max_ts;max_fs;max_K0;max_Q10];
diff_par=Max-Min;
J_old = 3000000;
for simu=1:nsimu1
		while (true)
            %propose new sets of parameters
        par_new = par_old+(rand(13,1)-0.5).*diff_par/allow;
        t = par_new(4)-par_new(5)*0.01*(100-CELLSAND_meanovertime);
        f31=par_new(6);
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
             && par_new(11)>Min(11)&&par_new(11)<Max(11)...
             && par_new(12)>Min(12)&&par_new(12)<Max(12)...
             && par_new(13)>Min(13)&&par_new(13)<Max(13)...
             && min(t(:))>0 && min(f21(:))>0)             
               break;  
            end
              
        end
         
        
   
%%%%%%%%%%%%%%% parameter space %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

D1=par_new(1);
D2=par_new(2);
zt=par_new(3);
fs=par_new(6:9);
K0=par_new(10:12);
Q10=par_new(13);
T_SCALAR=Q10.^((Tsoil-25)/10);
   for layer=1:10;  
        D_SCALAR(layer)=exp(-zsoi(layer)/zt);
    end;
    
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%% calculating a,b&c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % add advection later on
% calculate tridiagonal matrix for non-permofrost soil
diffus=D1;            
for s=1:3;
for j=1:10;   
    
    if j==1;
        d_m1_zm1(j)=0;
        w_p1(j)=(zsoi(j+1)-zisoi(j))/dz_node(j+1); % zsoi: node depth
        d_p1(j)=1/((1-w_p1(j))/diffus+w_p1(j)/diffus)/dz(j);  
        d_p1_zp1(j)=d_p1(j)/dz_node(j+1);
        f_m1(j)=0;%adv_flux(j);
        f_p1(j)=0;%adv_flux(j+1);
        pem1(j)=0;
        pep1(j)=f_p1(j)/d_p1_zp1(j);
    elseif j==10;
        w_m1(j)=(zisoi(j-1)-zsoi(j-1))/dz_node(j);
        d_m1(j)=1/((1-w_m1(j))/diffus+w_m1(j)/diffus)/dz(j);
        d_m1_zm1(j)=d_m1(j)/dz_node(j);
        d_p1_zp1(j)=0;
        f_m1(j)=0;%adv_flux(j);
        f_p1(j)=0;
        pem1(j)=f_m1(j)/d_m1_zm1(j);
%         pep1(j)=f_p1(j)/d_p1_zp1(j);
        pep1(j)=0;%f_p1(j)/d_p1_zp1(j); couldnot use the original one due to d_p1_zp1=0;
    else
        w_m1(j)=(zisoi(j-1)-zsoi(j-1))/dz_node(j);
        w_p1(j)=(zsoi(j+1)-zisoi(j))/dz_node(j+1);
        d_m1(j)=1/((1-w_m1(j))/diffus+w_m1(j)/diffus)/dz(j);
        d_p1(j)=1/((1-w_p1(j))/diffus+w_p1(j)/diffus)/dz(j);    
        f_m1(j)=0;%adv_flux(j);
        f_p1(j)=0;%adv_flux(j+1);
        d_m1_zm1(j)=d_m1(j)/dz_node(j);
        d_p1_zp1(j)=d_p1(j)/dz_node(j+1);
        pem1(j)=f_m1(j)/d_m1_zm1(j);
        pep1(j)=f_p1(j)/d_p1_zp1(j);
    end;
    
    aaa_pem1(j)=max(0,(1-0.1*abs(pem1(j)))^5);
    aaa_pep1(j)=max(0,(1-0.1*abs(pep1(j)))^5);    
         
   if j>0 && j<11;
    a_p_0(j)=0;
    a_tri(j)=-(d_m1_zm1(j)*aaa_pem1(j)+max(f_m1(j),0));
    c_tri(j)=-(d_p1_zp1(j)*aaa_pep1(j)+max(-f_p1(j),0));
    b_tri(j)=-a_tri(j)-c_tri(j)+a_p_0(j);
    record_a_tri((s-1)*10+j)=a_tri(j);
    record_c_tri((s-1)*10+j)=c_tri(j);
    record_b_tri((s-1)*10+j)=b_tri(j);
   end
end;
end
record_a_tri=record_a_tri';% left-off diagonal 
record_b_tri=record_b_tri';% diagonal
record_c_tri=record_c_tri';% right-off diagonal
record_abc_tri=[record_a_tri record_b_tri record_c_tri];
a=record_a_tri(2:30);b=record_b_tri;c=record_c_tri(1:29);
tridiag1 = gallery('tridiag',a,b,c);%tridiag1 is for non-permosfrost soil
% calculate tridiagonal matrix for permofrost soil
diffus=D2;            
for s=1:3;
for j=1:10;   
    
    if j==1;
        d_m1_zm1(j)=0;
        w_p1(j)=(zsoi(j+1)-zisoi(j))/dz_node(j+1); 
        d_p1(j)=1/((1-w_p1(j))/diffus+w_p1(j)/diffus)/dz(j);  
        d_p1_zp1(j)=d_p1(j)/dz_node(j+1);
        f_m1(j)=0;
        f_p1(j)=0;
        pem1(j)=0;
        pep1(j)=f_p1(j)/d_p1_zp1(j);
    elseif j==10;
        w_m1(j)=(zisoi(j-1)-zsoi(j-1))/dz_node(j);
        d_m1(j)=1/((1-w_m1(j))/diffus+w_m1(j)/diffus)/dz(j);
        d_m1_zm1(j)=d_m1(j)/dz_node(j);
        d_p1_zp1(j)=0;
        f_m1(j)=0;
        f_p1(j)=0;
        pem1(j)=f_m1(j)/d_m1_zm1(j);
        pep1(j)=0;
    else
        w_m1(j)=(zisoi(j-1)-zsoi(j-1))/dz_node(j);
        w_p1(j)=(zsoi(j+1)-zisoi(j))/dz_node(j+1);
        d_m1(j)=1/((1-w_m1(j))/diffus+w_m1(j)/diffus)/dz(j);
        d_p1(j)=1/((1-w_p1(j))/diffus+w_p1(j)/diffus)/dz(j);    
        f_m1(j)=0;
        f_p1(j)=0;
        d_m1_zm1(j)=d_m1(j)/dz_node(j);
        d_p1_zp1(j)=d_p1(j)/dz_node(j+1);
        pem1(j)=f_m1(j)/d_m1_zm1(j);
        pep1(j)=f_p1(j)/d_p1_zp1(j);
    end;
    
    aaa_pem1(j)=max(0,(1-0.1*abs(pem1(j)))^5);
    aaa_pep1(j)=max(0,(1-0.1*abs(pep1(j)))^5);    
         
   if j>0 & j<11;
    a_p_0(j)=0;
    a_tri(j)=-(d_m1_zm1(j)*aaa_pem1(j)+max(f_m1(j),0));
    c_tri(j)=-(d_p1_zp1(j)*aaa_pep1(j)+max(-f_p1(j),0));
    b_tri(j)=-a_tri(j)-c_tri(j)+a_p_0(j);
    record_a_tri((s-1)*10+j)=a_tri(j);
    record_c_tri((s-1)*10+j)=c_tri(j);
    record_b_tri((s-1)*10+j)=b_tri(j);
   end
end;
end
record_a_tri=record_a_tri';% left-off diagonal 
record_b_tri=record_b_tri';% diagonal
record_c_tri=record_c_tri';% right-off diagonal
record_abc_tri=[record_a_tri record_b_tri record_c_tri];
a=record_a_tri(2:30);b=record_b_tri;c=record_c_tri(1:29);
tridiag2 = gallery('tridiag',a,b,c);%tridiag2 is for permosfrost soil
%tridiag2=tridiag2';

%%%%%%%%%%%%%%% for equilibrium calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ph=0.5; % this is a place holder (ph) for f21 and will be replaced in the next step.
    aa=[-1 fs(2) fs(4);ph -1 0;fs(1) fs(3) -1]; % the transfer matrix
%%%%%%%%%%% A %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 soilcarbon_ss_emulator1=nan(288,192,30);
 soilcarbon_ss_emulator2=nan(288,192,30); 
k=1;
for lon=1:288;
    for lat=1:192;
                  %%%%%%%%%%% K %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            D_SCALAR_vr(lon,lat,1:10)=D_SCALAR';
            K_s1=K0(1)*(scalar_soil_mean(lon,lat,1:10)).*D_SCALAR_vr(lon,lat,1:10).*T_SCALAR(lon,lat,1:10); 
            K_s2=K0(2)*(scalar_soil_mean(lon,lat,1:10)).*D_SCALAR_vr(lon,lat,1:10).*T_SCALAR(lon,lat,1:10); 
            K_s3=K0(3)*(scalar_soil_mean(lon,lat,1:10)).*D_SCALAR_vr(lon,lat,1:10).*T_SCALAR(lon,lat,1:10); 
            K_s(lon,lat,1:10)=K_s1;
            K_s(lon,lat,11:20)=K_s2;
            K_s(lon,lat,21:30)=K_s3;
       if  soilcarbon_o_total_NCSCD(lon,lat)>0 && sum(I(lon,lat,1:20))>0 && prod(K_s(lon,lat,1:30))>0;
            k=k+1;
            A=-eye(30);
            for i=1:3;
                for j=1:3;
                    if aa(i,j)>0;
                    receivorP=i;donorP=j;
                        for depth=1:10;
                            aa(2,1)=f21(lon,lat,depth);
                            A((receivorP-1)*10+depth,(donorP-1)*10+depth)=aa(receivorP,donorP);
                        end;
                    end;
                end;
            end;       
            soilcarbon_ss_emulator1(lon,lat,1:30)=(A*diag(squeeze(K_s(lon,lat,:)))-tridiag1)\(squeeze(-I(lon,lat,:)));
            soilcarbon_ss_emulator2(lon,lat,1:30)=(A*diag(squeeze(K_s(lon,lat,:)))-tridiag2)\(squeeze(-I(lon,lat,:)));
        end
    end
end
soilcarbon_ss_emulator1(mask_perm2==1)=0;soilcarbon_ss_emulator2(mask_perm2==0)=0;
soilcarbon_ss_emulator=soilcarbon_ss_emulator1+soilcarbon_ss_emulator2;
soilcarbon_ss_emulator(mask_overall(:,:,:)==0)=NaN;
soilcarbon1_M=soilcarbon_ss_emulator(:,:,1:10).*soildepth_gridded(:,:,1:10);
soilcarbon2_M=soilcarbon_ss_emulator(:,:,11:20).*soildepth_gridded(:,:,1:10);
soilcarbon3_M=soilcarbon_ss_emulator(:,:,21:30).*soildepth_gridded(:,:,1:10);
soilcarbon123_M=soilcarbon1_M+soilcarbon2_M+soilcarbon3_M;
soilcarbon_M_top=sum(soilcarbon123_M(:,:,1:5),3)+soilcarbon123_M(:,:,6)*((0.3-zisoi(5))/dz(6));%0-30cm
soilcarbon_M_total=sum(soilcarbon123_M(:,:,1:7),3)+soilcarbon123_M(:,:,8)*((1-zisoi(7))/dz(8));%0-100cm
soilcarbon_M_100_200=soilcarbon123_M(:,:,8)*((zisoi(8)-1)/dz(8))+soilcarbon123_M(:,:,9)*((2-zisoi(8))/dz(9));% 100-200 cm
soilcarbon_M_200_300=soilcarbon123_M(:,:,9)*((zisoi(9)-2)/dz(9))+soilcarbon123_M(:,:,10)*((3-zisoi(9))/dz(10));% 100-200 cm
soilcarbon_M_top(soilcarbon_M_top==0)=NaN;
soilcarbon_M_total(soilcarbon_M_total==0)=NaN;
soilcarbon_M_100_200(soilcarbon_M_100_200==0)=NaN;
soilcarbon_M_200_300(soilcarbon_M_200_300==0)=NaN;
% logged scheme: changed
log_soilcarbon_M_top=log10(soilcarbon_M_top);
log_soilcarbon_M_total=log10(soilcarbon_M_total);
log_soilcarbon_M_100_200=log10(soilcarbon_M_100_200);
log_soilcarbon_M_200_300=log10(soilcarbon_M_200_300);
scalar_sd=0.5;
diff1=log_soilcarbon_M_top-log_soilcarbon_o_top_NCSCD;
J1=(diff1.^2)./(2*(scalar_sd*log_soilcarbon_o_top_NCSCD).^2);
diff2=log_soilcarbon_M_total-log_soilcarbon_o_total_NCSCD;
J2=(diff2.^2)./(2*(scalar_sd*log_soilcarbon_o_total_NCSCD).^2);
diff3=log_soilcarbon_M_100_200-log_soilcarbon_o_NCSCD200;
J3=(diff3.^2)./(2*(scalar_sd*log_soilcarbon_o_NCSCD200).^2);
diff4=log_soilcarbon_M_200_300-log_soilcarbon_o_NCSCD300;
J4=(diff4.^2)./(2*(scalar_sd*log_soilcarbon_o_NCSCD300).^2);
J_new=nansum(J2(:))+nansum(J3(:))+nansum(J4(:));
delta_J = J_new-J_old;
if min(1,exp(-delta_J))>rand
   Parameters_keep(:,upgrade)=par_new;
   J_keep(upgrade)=J_new;
   upgrade=upgrade+1
   par_old=par_new;
   J_old=J_new;   
   record_k(upgrade)=k;
end
simu
upgrade
Parameters_rec(:,simu)=par_old;      
J_rec(:,simu)=J_old;                                  
end