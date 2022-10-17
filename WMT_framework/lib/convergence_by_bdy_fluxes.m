function [drho_F] = convergence_by_bdy_fluxes(Qsm, Qsh, SWR, SA, rho, alpha, beta, dz, cp0)

% By Sjoerd Groskamp, April 2018.

% Here we calculate the CONVERGENCE of heat and freshwater/Salt in a grid-box, by  surface fluxes specifically.
% The input E,P,R (in m/s) and LW, SW, LH and SH (in W/m2).
% The T-grid is the middle of a box, where tracers are given.
%
% Qsm = Surface mass flux.
% Qsh = Surface heat flux.
% SWR = Shortwave radiation (3D).

%% Preperations
[xl,yl,zl,tl]       = size(SA);
load("Geothermal.mat", 'Geothermal'); % Load geothermal heat flux (DeLavergne 2016)

% Here we add one row to the SWR, so that its convergence 
% (taking a difference in z-direction) has the same dimensions as SA.
% This needs to be done, because SWR is given at the 
% grid-box interface, while the convergence is given at the T-grid.
% (like SA). We just add a zero at the bottom, so it does no harm, and it remains 
SWR1                = zeros(xl, yl, zl+1, tl); 
SWR1(:,:,1:end-1,:) = SWR;
SWR                 = SWR1; clear SWR1;

%% HEATFLUX: FROM: W/m2 TO rho*kg/s
% 1) - We multiply the heatflux (W/m2) with (alpha/cp0), to leave kg/m2/s.
% 2) - We take the convergence with depth. 
%      Meaning we get the value at the middle of the water column, AT the tracer-grid point.
%      For the air-sea flux, that means the flux enters the volume of the
%      first grid point, and will be divided by its depth.
%      For the SWR this is the case for each grid point.
%      Note, vol = dx * dy * dz = Axy * dz, so we here divide by dz, 
%      and later multiply by volume, instead of directly multiplying by just Axy.
%      This is for consistency between diffusive flux and surface flux, all
%      expressed in rho / s at the T-grid.
dz4d            = repmat(dz,[1,1,1,tl]);
Qfac            = -alpha ./ cp0 ./ dz4d; clear alpha
Qfac_surf       = squeeze(Qfac(:,:,1,:));

% TURN THE HEATFLUXES INTO MASS FLUXES: 
SHF_rho_kg_s	= Qsh .* Qfac_surf; % rho /s
SWR_rho_kg_s	= -diff(SWR ,[],3) .* Qfac;
BHF_rho_kg_s    = Geothermal .* Qfac; % rho/s

%% FRESHWATER FLUX: FROM: m/s to rho/s
% 1) - We multiply the freshwater flux (m/s) rho, to get kg/m2/s, as in the review.
% 2) - We take the convergence for each bin, which in this case, is just
%      the convergence over the surface bin. We do this, by dividing by dz.
% 3) - We multiply Qm with SA * beta, which is unitless, this still kg/m2/s.
SMF_rho_kg_s =  Qsm .* squeeze((rho(:,:,1,:) + 1000) .* beta(:,:,1,:) .* SA(:,:,1,:) ./ dz(:,:,1,:)); % in rho kg/s.

%% COMBINE HEAT AND FRESHWATER FLUX
drho_F            = SWR_rho_kg_s;
drho_F(:,:,1,:)   = squeeze(drho_F(:,:,1,:)) + SHF_rho_kg_s + SMF_rho_kg_s; % Air-Sea flux in rho/s on T-grid.

topo = zeros(size(drho_F));
topo(isnan(drho_F)) = 1;
[I] = find(cumsum(topo, 3) == 1) - xl * yl; % find indices of cells just above the sea floor
I = I(I>0); % ignore negative indices, they probably correspond to surface land
drho_F(I)   = drho_F(I) + BHF_rho_kg_s(I); % Sea-floor flux in rho/s on T-grid.