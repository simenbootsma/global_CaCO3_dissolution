function [SW_O03, SW_SF] = WMTR_SW_paramterisations(zt, SW, Chla, topo)

%% INPUT
% SW   = lon by lat by time, 3D matrix of Short Wave radiation at the surface, for different time steps.

%% Preallocation and stuff:
[xl,yl,tl]     = size(SW); % Sizes.
dim            = length(Chla);
if dim == 1
    disp(['CHLOROPHYL is constant!  The value is: ',num2str(Chla),' mg/m^3.'])
    Chla            = ones(xl,yl,tl) .* Chla; % Set Cholorphyl to minimum value, everywhere. Was 0.0158.
end

%% CALCULATE SW PENETRATION EFFECTS, REALISTIC
% The amount of heat converging in a certain part of the ocean, depends on the flux above and below it.
% The shortwave flux, at the surface, must always be 1.
% Tracers are given at inputs zt.
% We must fine the edges around zt, making sure that the first edge is the surface and is 1, the last edge is 0 is zero, and that the given depths zt, lies in between these surfaces.
% Note the bottom point is sacrifiesed
if zt(1) ~ 0; % If zt(1) is not at the surface:
   zt_edge = 0.5 * (zt(1:end-1) + zt(2:end)); % edges between zt. Length: zl-1;
   z       = [0,zt_edge]; % First edge must be zero, second edge must be below first zt. Length: zl.
elseif zt(1)==0 % If the first zt is zero.
   disp('DOUBLE CHECK: FIRST TRACER POINT IS AT SURFACE, THIS MAY GIVE PROBLEMS LATER ON')
   zt_edge = 0.5 * (zt(1:end-1) + zt(2:end)); % edges between zt-points, but only below the surface. Length: zl-1
   z       = [0,zt_edge(2:end)]; %First point at surface, but first tracer point, is now the first point below the surface, which                                   is the SECOND point of zt. LEngth: zl-1
end
zl                = length(z); % Number of depth steps.
z4d           	  = permute(repmat(repmat(z',[1,tl]),[1,1,xl,yl]),[3,4,1,2]); % 4D Depth values.
SW4d              = permute(repmat(SW,[1,1,1,zl]),[1,2,4,3]); % Repeat SW values, with depth.
SW4d              = SW4d + topo; % 4D SW values.
SW4d(isnan(SW4d)) = 0; % No SWR in the bottom.

% Making CHLA into 4D
% - Chla is cst with depth for MA94 (Page 1657 of their paper).
% - I ASSUME the same for O03, but I'm not 100% sure.
Chla4d = permute(repmat(Chla,[1,1,1,zl]),[1,2,4,3]); % 4D Chla, 

%% Ohlman 2003:
Chla4d(Chla4d>3)   = 3; % Valid for the range 0.01 - 3 mg/m^3, only 2% is larger, none is lower.
[A1,A2,B1,B2]      = O03_coefficents(Chla4d);
fac_O03            = A1 .* exp(-B1 .* z4d) + A2 .* exp(-B2 .* z4d);
idx1               = z4d<2;
fac_O03(idx1)      = 1; clear idx1 % The above is valid for depths grester then 2m, above that, it receives full stuff.
SW_O03             = SW4d .* fac_O03;
clear A1 A2 B1 B2 fac_O03

%% SURFACE ABSORPTION, CHECK CASE
fac_SF          = zeros(size(SW4d));
fac_SF(:,:,1,:) = ones(xl,yl,tl);
SW_SF           = SW4d .* fac_SF; clear SW4d