function [Rrange, Grho, Gtracer, dia_M] = run_framework(tracer)
% INPUT:
%   -> tracer: 4D (xyzt) matrix of dimensions 360x180x101x12 containing tracer
%   concentrations

% OUTPUT:
%   -> Grho: 2D matrix of dimensions 4xN containing the (bdy, N, H,
%   V)-components of the water mass transformation as a function of the
%   neutral density.
%   -> Gtracer: 2D matrix of dimensions 4xN containing the (bdy, N, H,
%   V)-components of the dia-boundary tracer transport as a function of the
%   neutral density.
%   -> dia_M: 2D matrix of dimensions 3xN containing the (N, H,
%   V)-components of the material tracer changes as a function of the
%   neutral density.

folder = "";
addpath("lib\");
tic

%% LOADING DATA
load(folder + "WOA18_gsw_N2LL_plus.mat"); % Load WOA18 WMT data
load(folder + "WOA18_gsw_N2LL_plus.mat", 'alpha', 'beta', 'gamma'); % Load WOA18 WMT data
addpath('gsw/');
addpath('gsw/library/');

% Clear Variables:
clear gamma_original ShallowMask Rh0 Rh0_surf SA_surf CT_surf rho_surf 
clear b_surf beta_surf alpha_surf Cb_surf Tb_surf
clear Rh0  N2_surf
clear LH_OA SH_OA LW_OA SW_OA E_OA P_OA R_OA 
clear albedo_diffusive albedo_direct albedo_seferian_noCHLA albedo_whc cloud_cover  
clear CTmax CTmin deg2m dir dR DR lat_A lon_A ml_pres_max ml_pres_min 
clear rho0 rho_fresh Rl Rmax Rmin Rrange SAmax SAmin savingtime Slope_max time_stamp 

%% Density Bins
% The density difference was set at a very small constant value in order to define and stabalise gamma.
% For the purpose here, we can actually use a larger number, and make it change.
Rmin        = min(gamma(:))-0.1; %Minimum Density (from WOA_STEP2.m) is 19
Rmax        = max(gamma(:))+0.1; %Maximum Density (from WOA_STEP2.m) is 28.3
dR          = 0.05;% Changes to 0.08, instead on 0.025.
Redge       = (Rmin-0.5*dR):dR:(Rmax+0.5*dR); % The edges of this range.
Rrange      = 0.5*(Redge(1:end-1)+Redge(2:end));
Rl          = length(Redge); % Length of the Transformation vector.
dSig        = diff(Redge);
wind        = 2; % window size for running average


%% Deal with mixed layer:
mlp_max                     = 1000;        % Set a maximum MLP (dbar).
mlp(mlp>mlp_max)            = mlp_max;


%% SELECT VARIABLES WE NEED IN THIS SCRIPT AND REMOVE THE REST TO %save MEMORY
gamma     	= gamma    	+ topo;
p         	= p         + topo;
b        	= b         + topo;
rho         = rho       + topo;
SA          = SA        + topo;
alpha      	= alpha    	+ topo;
beta        = beta      + topo;
tracer      = tracer    + topo;
gamma_surf  = gamma_surf + topo_surf;

% Surface fluxes:
LH_Core     = LH_Core   + topo_surf;
SH_Core     = SH_Core   + topo_surf;
LW_Core     = LW_Core   + topo_surf;
SW_Core     = SW_Core   + topo_surf;
E_Core      = E_Core    + topo_surf;
P_Core      = P_Core    + topo_surf;
R_Core      = R_Core    + topo_surf;


%% CREATE BALANCED PROFILE
% Surface mass and heat flux:
Qsm         = E_Core + P_Core + R_Core; % Surface Mass Flux.
Qsh         = LH_Core + SH_Core + LW_Core; % Surface Heat Flux.


%% SHORTWAVE PARAMTERIZATIONS
disp('Calculate SWR depth Penetration')
Xtopo               = zeros(size(topo));
Xtopo(:,:,1,:)      = topo_surf;
Xtopo(:,:,2:end,:)  = topo(:,:,1:end-1,:);
[SW_O03, ~]         = WMTR_SW_paramterisations(zt, SW_Core, Chla, Xtopo); clear Xtopo Chla


%% CONVERGENCE OF DENSITY DUE TO AIR-SEA FLUX (RHO / S):
dRho_F = convergence_by_bdy_fluxes(Qsm, Qsh, SW_O03, SA, rho, alpha, beta, dz, cp0); % rho / s on T-grid.


%% DEFINING DIFFUSIVITIES
[~, ~, zl, tl] = size(alpha);

% Constant diffusivities
% K                       = 500 .* ones(xl,yl,zl,tl); % Mesoscale Tracer diffusion in m2/s
% D                       = 1e-5 .* ones(xl,yl,zl,tl); % Small-scale Tracer diffusion in m2/s

% Variable diffusivitites
load('WOA18_K.mat', 'WOA18_K');
load('Casimir_D_Tgrid.mat', 'Casimir_D_Tgrid');
K = repmat(WOA18_K.K_WOA(:,:,2:end), [1, 1, 1, tl]);
D = Casimir_D_Tgrid.Dtot;


%% CONVERGENCE OF DENSITY DUE TO MIXING WITH CONSTANT K (RHO / S)
disp('CONVERGENCE WITH CONSTANT K')
load(folder + "WOA18_Slopes_interface_am.mat", 'Slopes_interface');
mlp4d                       = permute(repmat(mlp,[1,1,1,zl]),[1,2,4,3]);
Mdx                         = p<=mlp4d; % Mixed layer, where pressure is lower than mlp.
Ddx                      	= p>mlp4d;  clear mlp4d % Select where presure is larger than mlp (interior or deep ocean).
dRho_N                     	= convergence_by_isopycnal_mixing(alpha, beta, rho, K, dx, dy, dz, Slopes_interface); % rho / s
[dRho_H, dRho_V]           	= convergence_by_hor_and_ver_mixing(SA, CT, alpha, beta, rho, K, D, dx, dy, dz, zt); % rho / s
clear Slopes_interface;


%% BINNING FOR WMT (RHO/S, TO KG/S)
disp('BINNING WMT')
[Grho, Gtracer]   	= WMTR_binning_loop(dRho_N, dRho_H, dRho_V, dRho_F, tracer, gamma, rho, b, vol, Ddx, Mdx, Redge, wind);


%% TRACER FLUX AND BINNING:
disp('TRACER FLUX AND BINNING')
load(folder + "WOA18_Slopes_Tgrid_am.mat", 'Slopes_Tgrid'); 
[dTracer_N, dTracer_H, dTracer_V] = convergence_by_tracer_fluxes(tracer, rho, K, D, dx, dy, dz, vol, zt, Slopes_Tgrid); % g of tracer / m3 / s
clear Slopes_Tgrid;

dTracer_N = dTracer_N .* Ddx; % select interior
dTracer_H = dTracer_H .* Mdx; % select mixed layer

M_N       = binning_1D(dTracer_N,gamma,Redge) ./ tl;
M_H       = binning_1D(dTracer_H,gamma,Redge) ./ tl;
M_V       = binning_1D(dTracer_V,gamma,Redge) ./ tl;

M_N = RunningMean1D(M_N, wind);
M_H = RunningMean1D(M_H, wind);
M_V = RunningMean1D(M_V, wind);
dia_M = [M_N; M_H; M_V] ./ dSig;

toc;
end
