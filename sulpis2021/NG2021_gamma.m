%% Regional CaCO3 dissolution rates and bottom-up sinking flux estimates
%July 1st, 2020
%Code by Olivier Sulpis & Ashley Dinauer

%Please see associated paper:
%Sulpis O., Jeansson E., Dinauer A., Lauvset S. K., Middelburg, J. J. 
%2021
%Calcium carbonate dissolution patterns in the ocean
%Nature Geoscience


% ------------------------------------------------------------------------
% NOTE: 
% 
% This script is edited to work with neutral density instead of
% potential density, include negative values and use extended biomes that
% span the globe. For the original script and associated files, please go
% to: https://github.com/osulpis/pelagic_dissolution 
% ------------------------------------------------------------------------


clear all
close all

%% Load data files

tic
%load raw 2016 GLODAP dataset
G2_2021 = load('GLODAPv2_Merged_Master_File.mat');
%load raw 2016 GLODAP dataset (Olsen et al., 2016, ESSD)
%load seafloor CaCO3 burial rates (total sed burial rates are from Jahnke 1996 GBC, shared by Nicolas Gruber, and CaCO3 sed contents are from Jenkins 1997 ST)
%load Biomes from Fay and McKinley (2014 ESSD) and their 2D and 3D mapped versions
%load sinking CaCO3 flux sediment trap dataset (from Mouw et al., 2016, ESSD)
%load seawater age from TTD analysis (from Emil Jeansson)
%load seawater age from 14C (from Gebbie and Huybers 2012, JPO, shared by Adam Subhas)
load('starting_workspace1.mat')
load('starting_workspace2.mat')
load('starting_workspace3.mat')
load('starting_workspace4.mat')

% -----------------------------------------------------
% [Simen]
% -----------------------------------------------------
G2_2021.G2gamma(G2_2021.G2gamma == -9999) = nan;
I = find(~isnan(G2_2021.G2latitude) & ~isnan(G2_2021.G2longitude) & ~isnan(G2_2021.G2depth) & ~isnan(G2_2021.G2gamma));
F = scatteredInterpolant(G2_2021.G2latitude(I), G2_2021.G2longitude(I), G2_2021.G2depth(I), G2_2021.G2gamma(I), 'linear', 'none');
G2gamma0 = nan(size(G2latitude));
I2 = ~isnan(G2sigma0);
G2gamma0(I2) = F(G2latitude(I2), G2longitude(I2), G2depth(I2));
% -----------------------------------------------------

% -----------------------------------------------------
% [Simen]: extend biomes
% -----------------------------------------------------
I = isnan(G2biomes);
F2 = scatteredInterpolant(G2longitude(~I), G2latitude(~I), G2depth(~I), G2biomes(~I), 'nearest');
G2biomes(I) = F2(G2longitude(I), G2latitude(I), G2depth(I));
% -----------------------------------------------------

clear F F2;

%% Sort sediment-trap data by biome

%make the longitude of sediment trap data monotonically increasing
for i=1:size(M16lon)
    if M16lon(i)<0
       M16lon_monotonic(i,1)=M16lon(i)+360;
    else
        M16lon_monotonic(i,1)=M16lon(i);
    end
end

%assign a biome number to each sediment-trap data: interpolation to the nearest neighbour
M16biomes=interp2(lon,lat,squeeze(G2map_Biomes(:,:,1)),M16lon_monotonic,M16lat,'nearest');

%express the PIC flux measured by sediment traps in mol / m2 / a
M16picflux=M16picflux.*365.25.*1e-3./12;
M16picflux(M16duration<0)=NaN;

%% Compute seawater age

%make longitudes monotonic
for i=1:size(G2longitude)
    if G2longitude(i)<0
       G2longitude_monotonic(i,1)=G2longitude(i)+360;
    else
        G2longitude_monotonic(i,1)=G2longitude(i);
    end
end

%assign an 14C age to each GLODAP sample: linear interpolation
C14_age=interp3(GH_lat,GH_depth,GH_lon,GH_tmean,G2latitude,G2depth,G2longitude_monotonic,'linear');

%G2_age is the age that will be used in the analysis. 
% First, we state that the age of each GLODAP sample is that from the TTD analysis.
TTD_age=mageF12;
G2_age=TTD_age;
% G2_age(G2_age <= 30) = nan;

% Then we create a transition zone in which the used age is from a mixture of TTD and 14C
for i=1:size(TTD_age)
    if TTD_age(i) >= 200 && TTD_age(i) <= 300  %transition period between 200 and 300 years
        transition_coef(i) = (TTD_age(i)-200)./100;  %linear transition
        G2_age(i) = (1-transition_coef(i))*TTD_age(i) + transition_coef(i)*C14_age(i);
 elseif TTD_age(i) > 300 %Finally, above 300 years, we use 14C ages only 
        G2_age(i)=C14_age(i);
    end 
end

%% Compute Alk* excess alkalinity

%For samples with no measured TA but with measured pH and DIC, we compute TA
CO2SYSdata=CO2SYS(G2phtsinsitutp,G2tco2,3,2,G2salinity,G2temperature,G2temperature,G2pressure,G2pressure,G2silicate,G2phosphate,1,10,1);
G2talk_calc=CO2SYSdata(:,21)+2.*CO2SYSdata(:,22)+CO2SYSdata(:,24)+CO2SYSdata(:,25)+CO2SYSdata(:,26)+CO2SYSdata(:,27)-CO2SYSdata(:,28);
clear CO2SYSdata
G2talk_combined=NaN(size(G2talk,1),1);
G2talk_combined(~isnan(G2talk))=G2talk(~isnan(G2talk));
G2talk_combined(isnan(G2talk))=G2talk_calc(isnan(G2talk)); 

%Equations (1) and (5) from Carter et al., Biogeosciences, 2014
Alk_star=G2talk_combined+1.26.*G2nitrate-66.4.*G2salinity;

%TA*, computed from Feely et al. (2002)
TA_star= 0.5*(G2talk.*35./G2salinity - (148.7+(61.36*G2salinity)+(0.0941*(G2oxygen+170.*G2phosphate))-(0.582*G2theta)).*35./G2salinity) + 0.05922*G2aou;

%% Compute CaCO3 dissolution rates

%"1" to have bins centered around certain gamma0 values, "0" to have bins with constant nb of data pts
binningmethod = 1;

%find what the biome numbers are 
G2_b=unique(G2biomes(~isnan(G2biomes)))';

%Set up MonteCarlo
MC=2000; %number of Monte Carlo simulations
minimum_gamma_increment=0.002;
maximum_gamma_increment=0.05;
Alkstar_uncertainty=6; %[umol/kg], standard TAlk uncertainty from GLODAP
TAstar_uncertainty=6; %[umol/kg], standard TAlk uncertainty from GLODAP
age_uncertainty=0.2; %set to 20%, a standard relative uncertainty for CFC and 14C ages
minimum_gamma0=20;
maximum_gamma0=29;
max_gamma0bins=[minimum_gamma0:minimum_gamma_increment:maximum_gamma0]; 
%Results will be stored here: first dimension is Monte Carlo simulation number, 2nd dimension is biome number, 3rd dimension in density bin number 
Alkstar_rate=NaN(MC,size(G2_b,2),size(max_gamma0bins,2)); 
TAstar_rate=NaN(MC,size(G2_b,2),size(max_gamma0bins,2)); 
Alkstar_rate_depth=NaN(MC,size(G2_b,2),size(max_gamma0bins,2));
Alkstar_rate_maxdepth=NaN(MC,size(G2_b,2),size(max_gamma0bins,2));
Alkstar_rate_mindepth=NaN(MC,size(G2_b,2),size(max_gamma0bins,2));
Alkstar_rate_std_depth=NaN(MC,size(G2_b,2),size(max_gamma0bins,2));
Alkstar_rate_r2=NaN(MC,size(G2_b,2),size(max_gamma0bins,2));
Nbofpoints=NaN(MC,size(G2_b,2),size(max_gamma0bins,2));
Alkstar_range=NaN(MC,size(G2_b,2),size(max_gamma0bins,2));
Age_range=NaN(MC,size(G2_b,2),size(max_gamma0bins,2));

% -----------------------------------------------------
% [Simen]
% -----------------------------------------------------
Sigma_mean=NaN(MC,size(G2_b,2),size(max_gamma0bins,2)); 
Sigma_std=NaN(MC,size(G2_b,2),size(max_gamma0bins,2)); 
Gamma_mean=NaN(MC,size(G2_b,2),size(max_gamma0bins,2)); 
Gamma_std=NaN(MC,size(G2_b,2),size(max_gamma0bins,2)); 
% -----------------------------------------------------

for h=1:MC %FIRST FOR-LOOP: Monte Carlo for-loop 

    h %this is to display progression of for loop, can be commented out
    
    %Alk Star, age, and gamma increment will randomnly change at each MC simulation
    errAlkstar=2*Alkstar_uncertainty.*rand(size(G2_age,1),1);
    errTAstar=2*TAstar_uncertainty.*rand(size(G2_age,1),1);
    errage=G2_age.*2.*age_uncertainty.*rand(size(G2_age,1),1);
    gamma_increment=(maximum_gamma_increment-minimum_gamma_increment)*rand(1)+minimum_gamma_increment;
   %"random" values
    rAlk_star=Alk_star-Alkstar_uncertainty+errAlkstar;
    rTA_star=TA_star-TAstar_uncertainty+errTAstar;
    rG2_age=G2_age-G2_age.*age_uncertainty+errage;    
    rG2_age(rG2_age<0)=NaN;    
   
    for i=1:size(G2_b,2) %SECOND FOR-LOOP: for each biome, computes CaCO3 dissolution rates
    
        bgamma0=G2gamma0(G2biomes==G2_b(i));  %densities in a given biome
        bAlk_star=rAlk_star(G2biomes==G2_b(i));    %Alk_star in a given biome
        bTA_star=rTA_star(G2biomes==G2_b(i));      %TA_star in a given biome
        bage=rG2_age(G2biomes==G2_b(i));           %seawater age in a given biome
        bdepth=G2depth(G2biomes==G2_b(i));        %depths in a given biome
        blat=G2latitude(G2biomes==G2_b(i));        %lats in a given biome
        blon=G2longitude(G2biomes==G2_b(i));        %lons in a given biome
    
        [~,a]=sort(bgamma0,'descend'); %a is a dummy vector that sorts data according to descending gamma0
        bAlk_star=bAlk_star(a,:);    %Alk_star is sorted according to decreasing density
        bTA_star=bTA_star(a,:);    %Alk_star is sorted according to decreasing density
        bage=bage(a,:);                  %sw age is sorted according to decreasing density
        bgamma0=bgamma0(a,:);      %gamma0 is sorted according to decreasing density
        bdepth=bdepth(a,:);            %depth is sorted according to decreasing density
        blat=blat(a,:);            %lat is sorted according to decreasing density
        blon=blon(a,:);            %lon is sorted according to decreasing density

        % -----------------------------------------------------
        % [Simen]
        % -----------------------------------------------------
%         bgamma = G2gamma(G2biomes==G2_b(i));
%         bgamma = bgamma(a, :);
        % -----------------------------------------------------
        clear a
    
        a=~isnan(bAlk_star) & ~isnan(bTA_star) & ~isnan(bage) & ~isnan(bgamma0) & ~isnan(bdepth) & bdepth>100 & bage>0; 
        %a is a dummy vector that selects samples for which we have Alkstar, age, gamma, depth, positive age and below 100 m depth
        bAlk_star=bAlk_star(a==1);
        bTA_star=bTA_star(a==1);
        bage=bage(a==1);
        bgamma0=bgamma0(a==1);
        bdepth=bdepth(a==1);
        blat=blat(a==1);
        blon=blon(a==1);
        % -----------------------------------------------------
        % [Simen]
        % -----------------------------------------------------
%         bgamma = bgamma(a==1);
        % -----------------------------------------------------
    
        if binningmethod==0
        steps=30;  %steps sets how many data points we want in each density bin 
        else
        gamma0bins=[minimum_gamma0:gamma_increment:maximum_gamma0]; 
        end
    
        bin_nb=1;  %bin_nb is the number of the bin. We start at bin_nb = 1
        
        for j=1:size(gamma0bins,2) %THIRD FOR-LOOP
        %j=1:steps:size(bgamma0)-steps+1 %comment out if binning method=1
        
            if binningmethod==1
            bgamma0min=gamma0bins(j)-gamma_increment/2; %minimum gamma0 of the density bin
            bgamma0max=gamma0bins(j)+gamma_increment/2; %maximum gamma0 of the density bin
            bin_age=bage(bgamma0<=bgamma0max & bgamma0>=bgamma0min); %seawater ages within bin
            bin_Alk_star=bAlk_star(bgamma0<=bgamma0max & bgamma0>=bgamma0min); %Alkstar within bin
            bin_TA_star=bTA_star(bgamma0<=bgamma0max & bgamma0>=bgamma0min); %Alkstar within bin
            bin_gamma0=bgamma0(bgamma0<=bgamma0max & bgamma0>=bgamma0min); %denities within bin
%             bin_gamma = bgamma(bgamma0<=bgamma0max & bgamma0>=bgamma0min); % [Simen]
            bin_depth=bdepth(bgamma0<=bgamma0max & bgamma0>=bgamma0min); %depths within bin
            bin_lat=blat(bgamma0<=bgamma0max & bgamma0>=bgamma0min); %lat within bin
            bin_lon=blon(bgamma0<=bgamma0max & bgamma0>=bgamma0min); %lon within bin
            bin_nbsp=size(bin_depth,1); %number of data points within bin
            [b,~,~,~,stats]=regress(bin_Alk_star,[ones(size(bin_age)) bin_age]); %linear fit of ALk star versus age within
            [bb,~,~,~,~]=regress(bin_TA_star,[ones(size(bin_age)) bin_age]); %linear fit of TA star versus age within
            %bin. "b" is the slope of the fit, needs to be divided by two to
            %"bb" doesn't need to be divided be two because of the way TAstar is defined
            %give the dissolution rate. "stats" are stats associated with the fit
            else
            bin_age=bage(j:j+steps-1);                     %seawater age within the bin
            bin_Alk_star=bAlk_star(j:j+steps-1);       %Alk_star within the bin
            bin_TA_star=bTA_star(j:j+steps-1);       %Alk_star within the bin
            bin_gamma0=bgamma0(j:j+steps-1);         %gamma0 within the bin
%             bin_gamma=bgamma(j:j+steps-1);         %[Simen]
            bin_depth=bdepth(j:j+steps-1);               %depths within the bin
            bin_nbsp=size(bin_depth,1);
            %linear fit of Alk_star versus seawater age, intercept, slope and stats are saved
            [b,~,~,~,stats]=regress(bin_Alk_star,[ones(size(bin_age)) bin_age]);
            [bb,~,~,~,~]=regress(bin_TA_star,[ones(size(bin_age)) bin_age]);
            end
        
            %the following "if" is here to save only dissolution rates that fill some quality criteria, i.e., more than 10 data
            %points within the bin, age range within bin that is more than 20% of the mean age, Alk star range that is 
            %more than 12 umol/kg (which is twice the uncertainty)
            if bin_nbsp>10 && (max(bin_age)-min(bin_age))>0.2*nanmean(bin_age) && (max(bin_Alk_star)-min(bin_Alk_star))>12
            %binned_fits is the array where the results are saved


            % -----------------------------------------------------
            % [Simen]: plot fitting of this bin
            % -----------------------------------------------------
%             if bgamma0min > 28.4
%             figure; hold on
% %             plot(bin_age, bin_Alk_star, '.k', 'markersize', 20);
%             scatter(bin_age, bin_Alk_star, 20, bin_depth, 'filled')
% 
%             plot(bin_age, b(1) + b(2)*bin_age, '-r', 'linewidth', 2);
%             colormap(parula(55));
%             cb = colorbar;
%             caxis([0 5500])
%             cb.Label.String = 'Depth (m)';
%             cb.Label.FontSize = 14;
%             xlabel('Seawater age [years]')
%             ylabel('Alk* [umol/kg/year]')
%             title(string(round(bgamma0min, 4))+" \leq \sigma_0 \leq "+string(round(bgamma0max, 4))+", biome #"+string(i))
%             legend('data', "fit: rate = "+string(b(2))+" umol/kg/yr", 'location','best');
% 
%             show_world();hold on;
%             plot(bin_lon, bin_lat, '.r', 'markersize', 5);
% 
%             close all;
%             end
            % -----------------------------------------------------

            binned_fits(i,bin_nb,1)=b(1);                           %intercept of regression Alk_star versus age (umol / kg / a)
            binned_fits(i,bin_nb,2)=b(2);                           %slope of regression Alk_star versus age (umol / kg / a)
            binned_fits(i,bin_nb,3)=stats(1);                     %r squared of linear regression 
            binned_fits(i,bin_nb,4)=mean(bin_depth);      %mean depth of samples within a bin
            binned_fits(i,bin_nb,5)=std(bin_depth);          %standard deviation of depths within a bin
            binned_fits(i,bin_nb,6)=min(bin_depth);         %minimum depth within a bin
            binned_fits(i,bin_nb,7)=max(bin_depth);        %maximum depth within a bin
            binned_fits(i,bin_nb,8)=mean(bin_gamma0);   %mean gamma0 within a bin
            binned_fits(i,bin_nb,9)=std(bin_gamma0);       %standard deviation of gamma0 within a bin
            binned_fits(i,bin_nb,10)=min(bin_gamma0);    %minimum gamma0 within a bin
            binned_fits(i,bin_nb,11)=max(bin_gamma0);   %maximum gamma0 within a bin
            binned_fits(i,bin_nb,12)=bin_nbsp;                 %number of data points within a bin
            binned_fits(i,bin_nb,13)=mean(bin_Alk_star); %mean ALkstar 
            binned_fits(i,bin_nb,14)=std(bin_Alk_star);     %standard deviation of Alk star
            binned_fits(i,bin_nb,15)=min(bin_Alk_star);    %minimum Alk star
            binned_fits(i,bin_nb,16)=max(bin_Alk_star);   %maximum Alk star
            binned_fits(i,bin_nb,17)=mean(bin_age);        %mean age
            binned_fits(i,bin_nb,18)=std(bin_age);            %standard deviation of age
            binned_fits(i,bin_nb,19)=min(bin_age);           %minimum age
            binned_fits(i,bin_nb,20)=max(bin_age);          %maximum age 
            binned_fits(i,bin_nb,21)=bb(1);                       %intercept of regression Alk_star versus age (umol / kg / a)
            binned_fits(i,bin_nb,22)=bb(2);                       %slope of regression Alk_star versus age (umol / kg / a)
%             binned_fits(i,bin_nb,23)=mean(bin_gamma);   % [Simen]
%             binned_fits(i,bin_nb,24)=std(bin_gamma);    % [Simen]
%             binned_fits(i,bin_nb,25)=min(bin_gamma);    % [Simen]
%             binned_fits(i,bin_nb,26)=max(bin_gamma);    % [Simen]
            end
        
        bin_nb=bin_nb+1;    
        end
    
    end

    binned_fits(binned_fits==0)=NaN; 

    %important results are stored separately
    %1st dimension is Monte Carlo #, 2nd dimension is biome #, 3rd dimension is bin #
    Alkstar_rate(h,1:size(G2_b,2),1:size(binned_fits,2),:)=binned_fits(:,:,2)./2; %CaCO3 dissolution rate in [umol/kg/a]: slope of each fit divided by 2 (one mol of CaCO3 dissolved adds 2 moles of Alk)
    %/!\/!\ If TA_star is used instead of Alk_star, the slope is not divided by two
    TAstar_rate(h,1:size(G2_b,2),1:size(binned_fits,2),:)=binned_fits(:,:,22); %CaCO3 dissolution rate in [umol/kg/a]
    Alkstar_rate_depth(h,1:size(G2_b,2),1:size(binned_fits,2),:)=binned_fits(:,:,4);%mean depth for CaCO3 dissolution rate
    Alkstar_rate_maxdepth(h,1:size(G2_b,2),1:size(binned_fits,2),:)=binned_fits(:,:,7);%max depth for CaCO3 dissolution rate
    Alkstar_rate_mindepth(h,1:size(G2_b,2),1:size(binned_fits,2),:)=binned_fits(:,:,6);%min depth for CaCO3 dissolution rate
    Alkstar_rate_std_depth(h,1:size(G2_b,2),1:size(binned_fits,2),:)=binned_fits(:,:,5);%std dev of depth for CaCO3 dissolution rate
    Alkstar_rate_r2(h,1:size(G2_b,2),1:size(binned_fits,2),:)=binned_fits(:,:,3);%R2 of CaCO3 dissolution rate fit
    Nbofpoints(h,1:size(G2_b,2),1:size(binned_fits,2),:)=binned_fits(:,:,12);%nb of points for bin
    Alkstar_range(h,1:size(G2_b,2),1:size(binned_fits,2),:)=binned_fits(:,:,16)-binned_fits(:,:,15);%Alk star range for bin
    Age_range(h,1:size(G2_b,2),1:size(binned_fits,2),:)=binned_fits(:,:,20)-binned_fits(:,:,19);%Age range for bin
        
    % -----------------------------------------------------
    % [Simen]
    % -----------------------------------------------------
    Gamma_mean(h,1:size(G2_b,2),1:size(binned_fits,2),:)=binned_fits(:,:,8);
    Gamma_std(h,1:size(G2_b,2),1:size(binned_fits,2),:)=binned_fits(:,:,9);
    % -----------------------------------------------------

    gamma_range(h)=gamma_increment;%gamma increment used for a given Monte Carlo simulation
    clear binned_fits
%     Alkstar_rate(Alkstar_rate<0)=NaN;
%     TAstar_rate(TAstar_rate<0)=NaN;
    
end

save('Sulpis2021_data.mat', 'Alkstar_rate', 'Gamma_mean', 'Gamma_std');
disp("saved all data required for Simen's study!");


%% Compute seafloor accumulation and seafloor dissolution rate in continuous depth profiles

%a new continuous depth vector is created for interpolation: 6000 represents the maximum depth achievable
new_Alkstar_rate=NaN(MC,size(G2_b,2),6000);    %empty vector for interpolated dissolution rates 
new_TAstar_rate=NaN(MC,size(G2_b,2),6000);    %empty vector for interpolated dissolution rates 
new_Fpic=NaN(MC,size(G2_b,2),6000);   %empty vector for interpolated sinking fluxes (PIC fluxes)
new_acc=NaN(MC,size(G2_b,2),6000);     %empty vector for interpolated CaCO3 sediment accumulation rate at a given depth
new_swi=NaN(MC,size(G2_b,2),6000);      %empty vector for interpolated CaCO3 sediment dissolution rate at a given depth
new_wdepth=NaN(MC,size(G2_b,2),6000);%empty vector for interpolated depths
new_sfdepth=NaN(size(G2_b,2),6000); 
new_tdepth=[1:1:6000];

SWI_uncertainty=0.3.*SWI_Dcaco3;   %uncertainty on seafloor dissolution rate is set to 30% 
acc_uncertainty=0.3.*Acaco3;    %uncertainty on seafloor burial dissolution rate is set to 30%
sfdepth_uncertainty=200;    %uncertainty on seafloor depth is set to 200m

% FOR LOOP COMPUTING REGIONAL DEPTH PROFILE OF SEAFLOOR DISSOLUTION AND ACCUMULATION
for i=1:size(G2_b,2) 
    
            for h=1:MC

                dummysfdepth=Seafloor_depth;
                dummyacc=Acaco3;                    %CaCO3 seafloor burial rate [mol/m2/a]
                dummyswi=SWI_Dcaco3;           %CaCO3 seafloor dissolution rate [mol/m2/a]
                dummysfdepth=Seafloor_depth;  %seafloor depth

                dummysfdepth(G2map_Biomes~=i)=NaN;
                dummyacc(G2map_Biomes~=i)=NaN;
                dummyswi(G2map_Biomes~=i)=NaN;
                dummysfdepth(G2map_Biomes~=i)=NaN;
                
                randomstuff=rand(1);
                for randomlat=1:180
                    for randomlon=1:360
                         randommap(randomlat,randomlon)=randomstuff; %a map of random numbers each between 0 and 1
                    end
                end
                
                errSWI=2*SWI_uncertainty.*randommap;
                dummyswi=dummyswi-SWI_uncertainty+errSWI;    
                dummyswi(dummyswi<0)=1e-6;   %to avoid zeros in sediment dissolution rates, we set them to very small values
            
                erracc=2*acc_uncertainty.*randommap;
                dummyacc=dummyacc-acc_uncertainty+erracc;    
                dummyacc(dummyacc<0)=1e-6;   %to avoid zeros 
                          
                errdepth=2*sfdepth_uncertainty.*randommap;
                dummysfdepth=dummysfdepth-sfdepth_uncertainty+errdepth;    
                dummysfdepth(dummysfdepth<0)=1e-6;   %to avoid zeros 
                
                maxdepth=round(nanmax(nanmax((dummysfdepth))));
                mindepth=round(nanmin(nanmin((dummysfdepth))));
 
                if maxdepth>6000
                        maxdepth=6000;
                end
                if mindepth==0
                         mindepth=1;
                end
                new_sfdepth(i,301:maxdepth)=[301:1:maxdepth];
            
                baccumulation(i,:)=dummyacc(~isnan(dummyacc) & ~isnan(dummyswi) & ~isnan(dummysfdepth)); %CaCO3 burial rates for a given biome for every depth
                bswi(i,:)=dummyswi(~isnan(dummyacc) & ~isnan(dummyswi) & ~isnan(dummysfdepth));%CaCO3 sediment dissolution rates for a given biome for every depth
                bsfdepth(i,:)=dummysfdepth(~isnan(dummyacc) & ~isnan(dummyswi) & ~isnan(dummysfdepth)); 
     
                new_acc(h,i,:)=csaps(bsfdepth(i,:),baccumulation(i,:),1e-7,squeeze(new_sfdepth(i,:)));  %spline smoothing of CaCO3 burial rate
                new_swi(h,i,:)=csaps(bsfdepth(i,:),bswi(i,:),1e-7,squeeze(new_sfdepth(i,:)));  %spline smoothing of CaCO3 dissolution rate in sediments
            
                clear baccumulation bswi bsfdepth

            end
end

new_swi(new_swi<0)=1e-6;   %to avoid negative values
new_acc(new_acc<0)=1e-6;   %to avoid nagative values
new_sfsink=new_acc+new_swi; %sum of CaCO3 burial and sediment dissolution for a given biome

%% Compute water column dissolution rate and sinking rate in continuous depth profiles

%FIRST FOR-LOOP: biomes
for i=1:size(G2_b,2) %we select only the non-high latitudes biomes
    
    %SECOND FOR-LOOP: Monte Carlo simulations
    for h=1:MC
        
        if nansum(~isnan(Alkstar_rate(h,i,:)))>2 %we interpolate only if there are at least 2 dissolution rates in a given biome
          
            %Spline curve fitting, see https://nl.mathworks.com/help/curvefit/examples/cubic-smoothing-splines.html 
            %the nb of points in each bin is used as a weight for the spline curve fitting
            new_wdepth(h,i,round(min(Alkstar_rate_depth(h,i,:))):round(max(Alkstar_rate_depth(h,i,:))))=[round(min(Alkstar_rate_depth(h,i,:))):1:round(max(Alkstar_rate_depth(h,i,:)))];
            new_Alkstar_rate(h,i,:)=csaps(Alkstar_rate_depth(h,i,:),Alkstar_rate(h,i,:),5e-9,new_wdepth(h,i,:));%,Nbofpoints(h,i,:));    
            new_TAstar_rate(h,i,:)=csaps(Alkstar_rate_depth(h,i,:),TAstar_rate(h,i,:),5e-9,new_wdepth(h,i,:));%,Nbofpoints(h,i,:));     
            
            %compute sinking PIC flux at maximum depth for which a CaCO3
            %was derived in a given biome
            maxdepth=max(new_wdepth(h,i,:));
                                  
            SFsinkatmaxdepth=squeeze(new_sfsink(h,i,:));
            SFsinkatmaxdepth=SFsinkatmaxdepth(squeeze(new_sfdepth(i,:))==maxdepth);
            new_Fpic(h,i,maxdepth)=SFsinkatmaxdepth; %assign to that maximum depth a CaCO3 sediment accumulation rate
%             for j=maxdepth-1:-1:1 %integral from the maximum to the minimum depths of the CaCO3 dissolution rate times the vertical resolution
%                 new_Fpic(h,i,j)=new_rate(h,i,j)*1e-6*1029+new_Fpic(h,i,j+1); %in mol per m2 per a (1029 is a typical seawater density)               
%             end
            
        end
        
    end
    
end

%% Calcite and aragonite saturation depths

% dissolution above aragonite saturation depth / total dissolution
CO2SYSdata=CO2SYS(G2talk,G2phtsinsitutp,1,3,G2salinity,G2temperature,G2temperature,G2pressure,G2pressure,G2silicate,G2phosphate,1,10,1);
OmegaC=CO2SYSdata(:,30);
OmegaA=CO2SYSdata(:,31);
CO3=CO2SYSdata(:,22).*1e-6;
Ca=0.02128./40.087.*(G2salinity./1.80655); 
OmegaMgC=CO3.*Ca./(10^-5.89); %solubility of toadfish from Woosley et al. (2012, JGR) measured at T=25degC, atm. pressure, S=36.5 
clear CO2SYSdata

for i=1:size(G2_b,2)
    %calcite saturation depth: mean of depths with Omega C between 0.95 and 1.05
CSD(i,1)=nanmean(G2depth(G2biomes==i & OmegaC<1.05 & OmegaC>0.95));
CSDstd(i,1)=nanstd(G2depth(G2biomes==i & OmegaC<1.05 & OmegaC>0.95));

    %aragonite saturation depth: mean of depths with Omega A between 0.95 and 1.05
ASD(i,1)=nanmean(G2depth(G2biomes==i & OmegaA<1.05 & OmegaA>0.95));
ASDstd(i,1)=nanstd(G2depth(G2biomes==i & OmegaA<1.05 & OmegaA>0.95));

MgCSD(i,1)=nanmean(G2depth(G2biomes==i & OmegaMgC<1.05 & OmegaMgC>0.95));
MgCSDstd(i,1)=nanstd(G2depth(G2biomes==i & OmegaMgC<1.05 & OmegaMgC>0.95));

end

%% Regional plots
save('Sulpis2021_data_gamma_withNegatives_extendedBiomes_depth.mat', 'new_Alkstar_rate', 'new_tdepth');
toc
for i=1:size(G2_b,2) %we select only the non-high latitudes biomes
    
    figure(i)
    clf
    
    %plot CaCO3 dissolution rates with vertical error bars corresponding to the depth ranges
    subplot 121
    plot([0 0],[-6000 0],'k--','LineWidth',1.5)
    hold on
  
    %plot MC mean and 90% confidence intervals
    plot([-1,1],[-CSD(i)-CSDstd(i),-CSD(i)-CSDstd(i)],'Color',[0.4, 0.6, 1],'LineWidth',1)
    plot([-1,1],[-CSD(i)+CSDstd(i),-CSD(i)+CSDstd(i)],'Color',[0.4, 0.6, 1],'LineWidth',1)
    plut=plot([-1,1],[-CSD(i),-CSD(i)],'Color',[0.4, 0.6, 1],'LineWidth',2);
    
    plot([-1,1],[-ASD(i)-ASDstd(i),-ASD(i)-ASDstd(i)],'Color',[1, 0.6, 0.6],'LineWidth',1)
    plot([-1,1],[-ASD(i)+ASDstd(i),-ASD(i)+ASDstd(i)],'Color',[1, 0.6, 0.6],'LineWidth',1)
    plit=plot([-1,1],[-ASD(i),-ASD(i)],'Color',[1, 0.6, 0.6],'LineWidth',2);
    
    plot([-1,1],[-MgCSD(i)-MgCSDstd(i),-MgCSD(i)-MgCSDstd(i)],'Color',[255/255, 209/255, 25/255],'LineWidth',1)
    plot([-1,1],[-MgCSD(i)+MgCSDstd(i),-MgCSD(i)+MgCSDstd(i)],'Color',[255/255, 209/255, 25/255],'LineWidth',1)
    plut=plot([-1,1],[-MgCSD(i),-MgCSD(i)],'Color',[255/255, 209/255, 25/255],'LineWidth',2);
    
   plot1=plot(nanmean(squeeze(new_Alkstar_rate(:,i,:))),-new_tdepth,'Color',[0, 0, 0],'LineWidth',2);
   plot1=plot(nanmean(squeeze(new_Alkstar_rate(:,i,:)))-1.645.*nanstd(squeeze(new_Alkstar_rate(:,i,:))),-new_tdepth,'Color',[0, 0, 0],'LineWidth',1);
   plot1=plot(nanmean(squeeze(new_Alkstar_rate(:,i,:)))+1.645.*nanstd(squeeze(new_Alkstar_rate(:,i,:))),-new_tdepth,'Color',[0, 0, 0],'LineWidth',1);
%     plot1=plot(nanmean(squeeze(new_TAstar_rate(:,i,:))),-new_tdepth,'Color',[0, 0, 0],'LineWidth',2);
%     plot1=plot(nanmean(squeeze(new_TAstar_rate(:,i,:)))-1.645.*nanstd(squeeze(new_TAstar_rate(:,i,:))),-new_tdepth,'Color',[0, 0, 0],'LineWidth',1);
%     plot1=plot(nanmean(squeeze(new_TAstar_rate(:,i,:)))+1.645.*nanstd(squeeze(new_TAstar_rate(:,i,:))),-new_tdepth,'Color',[0, 0, 0],'LineWidth',1);
        
%     scatter(squeeze(Alkstar_rate(1,i,:)),-squeeze(Alkstar_rate_depth(1,i,:)),50,'filled','MarkerEdgeColor',[1 0.2 0.2],'MarkerFaceColor',[1 1 1])    
    xlabel('CaCO3 diss. r. (umol / kg /a)')
    ylabel('depth (m)')
    title(['Biome #',num2str(G2_b(i))])
    xlim([nanmin(nanmin(new_TAstar_rate(:,i,:))) nanmax(nanmax(new_TAstar_rate(:,i,:)))]) 
    ylim([-5500 0]) 
    xlim([-0.05 0.6]) 
    
    subplot 122
    
    hold on
 
    plot2=plot(nanmean(squeeze(new_Fpic(:,i,:))),-new_tdepth,'Color',[0, 0, 0],'LineWidth',2);
    plot2=plot(nanmean(squeeze(new_Fpic(:,i,:)))-1.645.*nanstd(squeeze(new_Fpic(:,i,:))),-new_tdepth,'Color',[0, 0, 0],'LineWidth',1);
    plot2=plot(nanmean(squeeze(new_Fpic(:,i,:)))+1.645.*nanstd(squeeze(new_Fpic(:,i,:))),-new_tdepth,'Color',[0, 0, 0],'LineWidth',1);
  
    %include sediment trap data
     bM16picflux=M16picflux(M16biomes==i);
     bM16depth=M16depth(M16biomes==i);
     bM16uniquedepth=unique(bM16depth);
     bM16duration=M16duration(M16biomes==i);
     
     %if several sediment traps are at the same depth and in the same biome, average them out 
     for k=1:size(bM16uniquedepth)
          bM16picmeanflux(k,1)=nanmean(bM16picflux(bM16depth==bM16uniquedepth(k)));
          bM16cumduration(k,1)=nansum(bM16duration(bM16depth==bM16uniquedepth(k)));
          bM16picflux_stddev(k,1)=nanstd(bM16picflux(bM16depth==bM16uniquedepth(k)));
     end
     e=errorbar(bM16picmeanflux,-bM16uniquedepth,bM16picmeanflux.*0,bM16picmeanflux.*0,bM16picflux_stddev,bM16picflux_stddev,'LineStyle','none');
     e.Color='black';
     e.CapSize=0;
     scatter(bM16picmeanflux,-bM16uniquedepth,'filled','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 0.2 0.2])
     clear bM16picmeanflux bM16picflux_stddev 

    xlim([0 1.3.*nanmax(nanmax(new_Fpic(:,i,:)))]) 
    ylim([-5500 0]) 
    title(['Biome #',num2str(G2_b(i))])
    xlabel('CaCO3 sinking flux (mol / m2 /a)')
    ylabel('depth (m)')
   xlim([0 0.8]) 
   
end
