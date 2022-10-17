function [dRho_H, dRho_V] = convergence_by_hor_and_ver_mixing(SA, CT, alpha, beta, rho, K, D, dx, dy, dz, zt)

% Sjoerd Groeskamp, April 2018.
% Here we calculate the convergence of heat and fresh water by 
% horiontal and vertical neutral diffusion.

%% HORIZONTAL AND VERTICAL GRADIENTS AT INTERFACE
[~,dTdx_Ayz,~,~,~,~,dTdy_Axz,~,~,~,~,dTdz_Axy,~,~,~] = gradients_interface_xyzt(CT,dx,dy,dz,zt);
[~,dSdx_Ayz,~,~,~,~,dSdy_Axz,~,~,~,~,dSdz_Axy,~,~,~] = gradients_interface_xyzt(SA,dx,dy,dz,zt);

 
%% Pre-calculations
[xl, yl, zl,tl]            = size(rho);
rho                         = rho + 1000;
[rho_Ayz, rho_Axz, rho_Axy] = tracer_on_interface_xyzt(rho);
[K_Ayz  ,K_Axz  ,~]         = tracer_on_interface_xyzt(K);
[~  ,~  ,D_Axy]             = tracer_on_interface_xyzt(D);

%% FLUXES AT INTERFACE

% In x-direction
F_CT_yz = rho_Ayz .* K_Ayz .* dTdx_Ayz; % Flux in CT * rho * m / s 
F_SA_yz = rho_Ayz .* K_Ayz .* dSdx_Ayz; % Flux in SA * rho * m / s  

% In y-direction
F_CT_xz = rho_Axz .* K_Axz .* dTdy_Axz;
F_SA_xz = rho_Axz .* K_Axz .* dSdy_Axz;

% In z-direction
F_CT_xy = rho_Axy .* D_Axy .* dTdz_Axy;
F_SA_xy = rho_Axy .* D_Axy .* dSdz_Axy;

%% FLUX DIVERGENCE AT T-GRID

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
dRho_H      = alpha .* (dFCT_dx + dFCT_dy) - beta .* (dFSA_dx + dFSA_dy); % Total, rho / s
dRho_V      = alpha .* dFCT_dz - beta .* dFSA_dz; % Total, rho / s