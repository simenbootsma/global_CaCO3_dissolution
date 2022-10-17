function [Grho, Gtracer] = WMTR_binning_loop(dRho_N, dRho_H, dRho_V, dRho_F, tracer, gamma, rho, b, vol, Ddx, Mdx, Redge, wind)

% Sjoerd Groeskamp, April 2018. [modified by Simen Bootsma, November 2021]
% Here we bin the density changes in density-bins to calculate WMT.
% The input is the net convergence of density, at the T-grid point.

%% Preperations
[~,~,~,tl] = size(gamma);

%% FLUXES: FROM: rho/s to rho kg /s
Fbdy    =        (rho + 1000) .* b .* dRho_F  .* repmat(vol,[1,1,1,tl]);
J_N   	= Ddx .* (rho + 1000) .* b .* dRho_N  .* repmat(vol,[1,1,1,tl]); % dRho_N (rho/s) * vol (m^3) * b (unitless) * rho = rho kg / s  
J_H   	= Mdx .* (rho + 1000) .* b .* dRho_H  .* repmat(vol,[1,1,1,tl]); % dRho_N (rho/s) * vol (m^3) * b (unitless) * rho = rho kg / s  
J_V   	=        (rho + 1000) .* b .* dRho_V  .* repmat(vol,[1,1,1,tl]); % dRho_N (rho/s) * vol (m^3) * b (unitless) * rho = rho kg / s  

%% ADDING Tracer
Tfac     	= tracer ./ (rho + 1000);      % g/m3 -> g/kg.
T_Fbdy      = Fbdy .* Tfac; % grams of Tracer * rho / s, due to air-sea fluxes.
T_J_N       = J_N  .* Tfac; % grams of Tracer * rho / s, due to Neutral Diffusion.
T_J_H       = J_H  .* Tfac; % grams of Tracer * rho / s, due to Horizontal Diffusion.
T_J_V       = J_V  .* Tfac; % grams of Tracer * rho / s, due to Vertical Diffusion.

% save('alk_buoyancyFluxes_varK_100mNaN.mat', 'T_Fbdy', 'T_J_N', 'T_J_H', 'T_J_V');

%% ACCUMARRAY MAY NOT DEAL WITH NAN
% Determine Indeces
dSig      	= diff(Redge);

%% BINNING
% BINNING THE WMT for density, Total:
G_bdy       = binning_1D(Fbdy,gamma,Redge) ./dSig / tl;
G_N         = binning_1D(J_N,gamma,Redge)  ./dSig / tl;
G_H         = binning_1D(J_H,gamma,Redge)  ./dSig / tl;
G_V         = binning_1D(J_V,gamma,Redge)  ./dSig / tl;

% BINNING THE WMT for density with Tracer, Total:
Gt_bdy       = binning_1D(T_Fbdy,gamma,Redge) ./dSig / tl;
Gt_N         = binning_1D(T_J_N,gamma,Redge)  ./dSig / tl;
Gt_H         = binning_1D(T_J_H,gamma,Redge)  ./dSig / tl;
Gt_V         = binning_1D(T_J_V,gamma,Redge)  ./dSig / tl;

%% Smooth the WMT

% WMT, Density, Total
G_bdy	= RunningMean1D(G_bdy,wind); % kilo gram (kg)
G_N    	= RunningMean1D(G_N,wind)    ; 
G_H  	= RunningMean1D(G_H,wind)    ; 
G_V    	= RunningMean1D(G_V,wind)  ; 

% WMT, with Tracer, Total
Gt_bdy	= RunningMean1D(Gt_bdy,wind); % (g)
Gt_N  	= RunningMean1D(Gt_N,wind)    ; 
Gt_H  	= RunningMean1D(Gt_H,wind)    ;      
Gt_V   	= RunningMean1D(Gt_V,wind)  ;  

% WMT totals:
Grho        = [G_bdy    ;G_N    ;G_H    ;G_V    ];
Gtracer     = [Gt_bdy   ;Gt_N   ;Gt_H   ;Gt_V   ];