function [dRho_N] = convergence_by_isopycnal_mixing(alpha, beta, rho, K, dx, dy, dz, Slopes_interface)

% Sjoerd Groeskamp, April 2018.
% Here we calculate the convergence of heat and fresh water by neutral
% diffusion.

%% Load Slopes
% Load the S, T and/or P gradients along neutral surfaces.
% The work here is based on 
% Groeskamp, Barker, McDougall, Abernathey, Griffies, that will be submitted to JPO in May 2018.
% These gradients are super accurate, more crude and inaccurate version 
% can also be used (but be warned of side effects described in the paper above).

dSAdx_N_Ayz = Slopes_interface.dSAdx_N_Ayz; % Salinity gradient, along neutral surface, in x-direction, AT the interface (!) I call Ayz (so dy * dz area)
dSAdy_N_Axz = Slopes_interface.dSAdy_N_Axz; % Salinity gradient, along neutral surface, in y-direction, AT the interface (!) I call Axz (so dx * dz area)
dSAdz_N_Axy = Slopes_interface.dSAdz_N_Axy; % Salinity gradient, along neutral surface, in z-direction, AT the interface (!) I call Axy (so dx * dy area)
dCTdx_N_Ayz = Slopes_interface.dCTdx_N_Ayz; % Temperature gradient, along neutral surface, in x-direction, AT the interface (!) I call Ayz (so dy * dz area)
dCTdy_N_Axz = Slopes_interface.dCTdy_N_Axz;
dCTdz_N_Axy = Slopes_interface.dCTdz_N_Axy;
clear Slopes_interface

%% Pre-calculations
[xl, yl, zl, tl]            = size(rho);
rho                         = rho + 1000;
[rho_Ayz,rho_Axz,rho_Axy]   = tracer_on_interface_xyzt(rho); % Tracers values on T-grid interfaces by averaging.
[K_Ayz  ,K_Axz  ,K_Axy]    	= tracer_on_interface_xyzt(K); % Tracers values on T-grid interfaces by averaging.

%% NEUTRAL FLUX AT INTERFACE

% In x-direction
F_CT_yz = rho_Ayz .* K_Ayz .* dCTdx_N_Ayz; % Flux in CT * rho * m / s  
F_SA_yz = rho_Ayz .* K_Ayz .* dSAdx_N_Ayz; % Flux in SA * rho * m / s

% In y-direction
F_CT_xz = rho_Axz .* K_Axz .* dCTdy_N_Axz;
F_SA_xz = rho_Axz .* K_Axz .* dSAdy_N_Axz;

% In z-direction
F_CT_xy = rho_Axy .* K_Axy .* dCTdz_N_Axy;
F_SA_xy = rho_Axy .* K_Axy .* dSAdz_N_Axy;

%% NEUTRAL FLUX DIVERGENCE AT T-GRID

% DIVERGENCE OF NEUTRAL FLUX, X-DIRECTION, AT T-GRID
dFCT_dx                    = (F_CT_yz - F_CT_yz([end,1:end-1],:,:,:))./repmat(dx,[1,1,1,tl]); % Divergence of Flux, CT * rho /s 
dFSA_dx                    = (F_SA_yz - F_SA_yz([end,1:end-1],:,:,:))./repmat(dx,[1,1,1,tl]); % Divergence of Flux, SA * rho /s 

% DIVERGENCE OF NEUTRAL FLUX, Y-DIRECTION, AT T-GRID
FCT_y_inbetween                 = NaN(xl,yl+1,zl,tl); % One extra, to calculate difference.
FCT_y_inbetween(:,2:end,:,:)    = F_CT_xz;
dFCT_dy                         = (FCT_y_inbetween(:,2:end,:,:) - FCT_y_inbetween(:,1:end-1,:,:))./repmat(dy,[1,1,1,tl]); clear FCT_y_inbetween

FSA_y_inbetween                 = NaN(xl,yl+1,zl,tl); % One extra, to calculate difference.
FSA_y_inbetween(:,2:end,:,:)    = F_SA_xz;
dFSA_dy                         = (FSA_y_inbetween(:,2:end,:,:) - FSA_y_inbetween(:,1:end-1,:,:))./repmat(dy,[1,1,1,tl]); clear FSA_y_inbetween

% DIVERGENCE OF NEUTRAL FLUX, Z-DIRECTION, AT T-GRID
FCT_z_inbetween                 = NaN(xl,yl,zl+1,tl); % One extra, to calculate difference.
FCT_z_inbetween(:,:,2:end,:)    = F_CT_xy;
dFCT_dz                         = (FCT_z_inbetween(:,:,1:end-1,:) - FCT_z_inbetween(:,:,2:end,:))./repmat(dz,[1,1,1,tl]); clear FCT_z_inbetween

FSA_z_inbetween                 = NaN(xl,yl,zl+1,tl); % One extra, to calculate difference.
FSA_z_inbetween(:,:,2:end,:)    = F_SA_xy;
dFSA_dz                         = (FSA_z_inbetween(:,:,1:end-1,:) - FSA_z_inbetween(:,:,2:end,:))./repmat(dz,[1,1,1,tl]); clear FSA_z_inbetween

%% DENSITY CHANGE AT T-GRID
dRho_N      = alpha .* (dFCT_dx + dFCT_dy + dFCT_dz) - beta .* (dFSA_dx + dFSA_dy + dFSA_dz); % Total, rho / s
