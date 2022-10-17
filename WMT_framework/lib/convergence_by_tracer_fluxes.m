function [dTracer_N, dTracer_H, dTracer_V] = convergence_by_tracer_fluxes(tracer, rho, K, D, dx, dy, dz, vol, zt, Slopes)
% Sjoerd Groeskamp, April 2018. [Modified by Simen Bootsma, November 2021]
% Here we calculate the convergence of some tracer by diffusion in 
% the neutral, horizontal and vertical direction.

%% Load Slopes
% Load the exact Neutral Slopes, based on the work of 
% Groeskamp, Barker, McDougall, Abernathey, Griffies, that will be submitted to JPO in May 2018.
% These slopes are super accurate, more crude and inaccurate version 
% can also be used (but be warned of side effects described in the paper above).
[Sx_Ayz, ~, Sx_Axy] = tracer_on_interface_xyzt(Slopes.Sx);
[~, Sy_Axz, Sy_Axy] = tracer_on_interface_xyzt(Slopes.Sy);
clear Slopes

%% Pre-calculations
[xl, yl, zl, tl]            = size(rho);
rho                         = rho + 1000;
[rho_Ayz,rho_Axz,rho_Axy]   = tracer_on_interface_xyzt(rho);
[K_Ayz  ,K_Axz  ,K_Axy]    	= tracer_on_interface_xyzt(K);
[~  ,~  ,D_Axy]             = tracer_on_interface_xyzt(D);
% tracer                      = tracer ./ (rho + 1000); % Tracer in g/m3 to g/kg, also same as saying Tracer/rho.
tracer                      = tracer ./ rho; % Tracer in g/m3 to g/kg, also same as saying Tracer/rho.

%% ISOPYCNAL GRADIENTS FROM SLOPES
% Tracer Gradients on interface:
[~,dTracerdx_Ayz,~,dTracerdx_Axy,...
 ~,~,dTracerdy_Axz,dTracerdy_Axy,...
 ~,dTracerdz_Ayz,dTracerdz_Axz,dTracerdz_Axy,~,~,~] = gradients_interface_xyzt(tracer,dx,dy,dz,zt);
   
% Neutral Gradients on interface:
S2          = Sx_Axy.^2 + Sy_Axy.^2;
dTracerdx_N_Ayz = dTracerdx_Ayz                                 + dTracerdz_Ayz .* Sx_Ayz; % 
dTracerdy_N_Axz = dTracerdy_Axz                                 + dTracerdz_Axz .* Sy_Axz;
dTracerdz_N_Axy = dTracerdx_Axy .* Sx_Axy + dTracerdy_Axy .* Sy_Axy + dTracerdz_Axy .* S2;
% clear Sx_Ayz Sx_Ayz Sy_Axz Sx_Axy Sy_Axy

% %% NEUTRAL FLUX AT INTERFACE
% F_Tracer_yz = rho_Ayz .* K_Ayz .* dTracerdx_N_Ayz;% In x-direction, Flux in rho * m2 / s * Tracer/rho/m = Tracer * m/s
% F_Tracer_xz = rho_Axz .* K_Axz .* dTracerdy_N_Axz;% In y-direction
% F_Tracer_xy = rho_Axy .* K_Axy .* dTracerdz_N_Axy;% In z-direction
% 
% %% NEUTRAL FLUX DIVERGENCE AT T-GRID
% dFTracer_dx_N = (F_Tracer_yz - F_Tracer_yz([end,1:end-1],:,:,:))./repmat(dx,[1,1,1,tl]); % Divergence of Flux, Tracer/s
% dFTracer_dy_N = (F_Tracer_xz - F_Tracer_xz(:,[end, 1:end-1],:,:))./repmat(dy,[1,1,1,tl]); % Divergence of Flux, Tracer/s
% dFTracer_dz_N = (F_Tracer_xy(:,:,[end,1:end-1],:) - F_Tracer_xy)./repmat(dz,[1,1,1,tl]); % Divergence of Flux, Tracer/s
% dFTracer_dy_N(:,1,:,:) = nan; % undefined (90N minus 90S)
% dFTracer_dz_N(:,:,1,:) = nan; % undefined (bottom minus surface)
% 
% % DIVERGENCE OF NEUTRAL FLUX, X-DIRECTION, AT T-GRID
% % dFTracer_dx_N                       = (F_Tracer_yz - F_Tracer_yz([end,1:end-1],:,:,:))./repmat(dx,[1,1,1,tl]); % Divergence of Flux, Tracer/s
% 
% % DIVERGENCE OF NEUTRAL FLUX, Y-DIRECTION, AT T-GRID
% % FTracer_y_inbetween                 = NaN(xl,yl+1,zl,tl); % One extra, to calculate difference.
% % FTracer_y_inbetween(:,2:end,:,:)    = F_Tracer_xz;
% % dFTracer_dy_N                       = (FTracer_y_inbetween(:,2:end,:,:) - FTracer_y_inbetween(:,1:end-1,:,:))./repmat(dy,[1,1,1,tl]); clear FTracer_y_inbetween
% 
% % DIVERGENCE OF NEUTRAL FLUX, Z-DIRECTION, AT T-GRID
% % FTracer_z_inbetween                 = NaN(xl,yl,zl+1,tl); % One extra, to calculate difference.
% % FTracer_z_inbetween(:,:,2:end,:)    = F_Tracer_xy;
% % dFTracer_dz_N                       = (FTracer_z_inbetween(:,:,1:end-1,:) - FTracer_z_inbetween(:,:,2:end,:))./repmat(dz,[1,1,1,tl]); clear FTracer_z_inbetween
% 
% clear F_Tracer_yz F_Tracer_xz F_Tracer_xy
% 
% %% HORIZONTAL AND VERTICAL FLUXES AT INTERFACE
% F_Tracer_yz = rho_Ayz .* K_Ayz .* dTracerdx_Ayz;% In x-direction
% F_Tracer_xz = rho_Axz .* K_Axz .* dTracerdy_Axz;% In y-direction
% F_Tracer_xy = rho_Axy .* D_Axy .* dTracerdz_Axy;% In z-direction
% 
% %% FLUX DIVERGENCE AT T-GRID
% % DIVERGENCE OF HORIZONTAL FLUX, X-DIRECTION, AT T-GRID
% dFTracer_dx                    = (F_Tracer_yz - F_Tracer_yz([end,1:end-1],:,:,:))./repmat(dx,[1,1,1,tl]); % Divergence of Flux, CT * rho /s 
% 
% % DIVERGENCE OF HORIZONTAL FLUX, Y-DIRECTION, AT T-GRID
% FTracer_y_inbetween                 = NaN(xl,yl+1,zl,tl); % One extra, to calculate difference.
% FTracer_y_inbetween(:,2:end,:,:)    = F_Tracer_xz;
% dFTracer_dy                         = (FTracer_y_inbetween(:,2:end,:,:) - FTracer_y_inbetween(:,1:end-1,:,:))./repmat(dy,[1,1,1,tl]); clear FTracer_y_inbetween
% 
% % DIVERGENCE OF VERTICAL FLUX, Z-DIRECTION, AT T-GRID
% FTracer_z_inbetween                 = NaN(xl,yl,zl+1,tl); % One extra, to calculate difference.
% FTracer_z_inbetween(:,:,2:end,:)    = F_Tracer_xy;
% dFTracer_dz                         = (FTracer_z_inbetween(:,:,1:end-1,:) - FTracer_z_inbetween(:,:,2:end,:))./repmat(dz,[1,1,1,tl]); clear FTracer_z_inbetween
% 
% %% Tracertrate Convergence on T-grid in grams of Tracertrate per second
% dTracer_N      = repmat(vol,[1,1,1,tl]) .* (dFTracer_dx_N + dFTracer_dy_N + dFTracer_dz_N); % Total, m3 * (Tracer / s) = m3 * (g of Tracer / m3 / s) = g of Tracer / s
% dTracer_H      = repmat(vol,[1,1,1,tl]) .* (dFTracer_dx   + dFTracer_dy);
% dTracer_V      = repmat(vol,[1,1,1,tl]) .* (dFTracer_dz);
% 

%% Alternative calculation [Simen 20-12-2021]
load('WOA18_gsw_N2LL_plus.mat', 'gamma');
[~,dgdx_Ayz,~,dgdx_Axy,...
 ~,~,dgdy_Axz,dgdy_Axy,...
 ~,dgdz_Ayz,dgdz_Axz,dgdz_Axy,~,~,~] = gradients_interface_xyzt(gamma,dx,dy,dz,zt);

% Neutral direction
dgdx_N_Ayz = dgdx_Ayz                                 + dgdz_Ayz .* Sx_Ayz;
dgdy_N_Axz = dgdy_Axz                                 + dgdz_Axz .* Sy_Axz;
dgdz_N_Axy = dgdx_Axy .* Sx_Axy + dgdy_Axy .* Sy_Axy + dgdz_Axy .* S2;

F_Tracer_N_x = rho_Ayz .* K_Ayz .* dTracerdx_N_Ayz .* dgdx_N_Ayz;  % rho * m^2/s * Tracer/kg/m * rho/m = rho * Tracer/m^3/s
F_Tracer_N_y = rho_Axz .* K_Axz .* dTracerdy_N_Axz .* dgdy_N_Axz;
F_Tracer_N_z = rho_Axy .* K_Axy .* dTracerdz_N_Axy .* dgdz_N_Axy;

% Horizontal and vertical directions
F_Tracer_x = rho_Ayz .* K_Ayz .* dTracerdx_Ayz .* dgdx_Ayz;
F_Tracer_y = rho_Axz .* K_Axz .* dTracerdy_Axz .* dgdy_Axz;
F_Tracer_z = rho_Axy .* D_Axy .* dTracerdz_Axy .* dgdz_Axy;

dTracer_N      = repmat(vol,[1,1,1,tl]) .* (F_Tracer_N_x + F_Tracer_N_y + F_Tracer_N_z); % Total, m3 * (rho * Tracer/m^3/ss) = rho * Tracer / s
dTracer_H      = repmat(vol,[1,1,1,tl]) .* (F_Tracer_x   + F_Tracer_y);
dTracer_V      = repmat(vol,[1,1,1,tl]) .* (F_Tracer_z);



