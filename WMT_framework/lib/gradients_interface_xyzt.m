function [dTdx,dTdx_Ayz,dTdx_Axz,dTdx_Axy,...
          dTdy,dTdy_Ayz,dTdy_Axz,dTdy_Axy,...
          dTdz,dTdz_Ayz,dTdz_Axz,dTdz_Axy,...
          T_Ayz,T_Axz,T_Axy] = gradients_interface_xyzt(tracer,dx,dy,dz,zt)
% This script calcultes the gradient of a tracer on interfaces that suround
% the T-grid or tracer grid.
% It is made for [x,y,z,t]-coordinates.

% Requirements and assumptions:
% - It is based on having Tracer values at the surface. zt(1) = 0m.

%% Pre calculations and explanations:
% Length of vectors.
[xl,yl,zl,tl]   = size(tracer);
% We do not want to devide by zero when calculating gradients.              
dx(dx==0)       = NaN; 
dy(dy==0)       = NaN;
dz(dz==0)       = NaN;

%% Explanation: the interfaces
% There are 3 gradients DCdx, dCdy and dCdz.
% There are 3 surfaces: - dy*dz, in x-direction = Ayz.
%                       - dx*dz, in x-direction = Axz
%                       - dx*dy, in x-direction = Axy
% Hence we have 9 results, the gradients in 3 directions on 3 surfaces.
%
% 1) We will calculate the tracer gradients at each tracer-grid point,
%       using the surrounding tracer values.
% 2) We will average these gradient-values on the surrounding surfaces
%       through which diffusion takes place.
%
% We first define the surface through which diffusion takes place:
%
% In x-direction we have Ayz interfaces.
% Because the world is 'round', we find xt = 0 = 360.
% Hence we have an interface, EAST of each T-grid.
% There are [xl,yl,zl] Ayz-interfaces.
%
% In y-direction we have Axz interfaces.
% We have an interface NORTH of each T-grid.
% The point at 90N/90S is a point. There is no diffusion through a point.
% We use: yt = 90S = []   (On antartica)
% We use: yt = 90N = NaN. (Yes, this is an assumption for simplicity)
% There are [xl,yl,zl] Axz-interfaces.
% 
% In z-direction we have Axy-interfaces.
% We have an interface ABOVE each T-grid.
% No diffusion through the surface or below the deepest T-grid (or bottom).
% We use that Axy(1) = at 0m = NaN.
% There are [xl,yl,zl] Axz-interfaces.
%
[T_Ayz,T_Axz,T_Axy] = tracer_on_interface_xyzt(tracer);
%% dCdx on each interface
% Average value of tracer, on interface between two tracer-values.
% Interface is to the East, in x-direction, on Ayz-interface:
%%%% T_Ayz                   = 0.5*(tracer + tracer([2:end,1],:,:,:));

% Resulting Tracer Gradient on T-grid.
dTdx                    = (T_Ayz-T_Ayz([end,1:end-1],:,:,:))./repmat(dx,[1,1,1,tl]); clear Tx_inbetween

% x-Gradient on Ayz Interface, in x-direction:
% Us that dx does not change, so we can use dx.
dTdx_Ayz                = (tracer([2:end,1],:,:,:)-tracer)./repmat(dx,[1,1,1,tl]); clear dx; 

% x-Gradient on Axz Interface:
dTdx_Axz                = NaN(xl,yl,zl,tl);
dTdx_Axz(:,1:end-1,:,:) = 0.5*(dTdx(:,1:end-1,:,:) + dTdx(:,2:end,:,:));

% x-Gradient on Axy Interface:
dTdx_Axy               	= NaN(xl,yl,zl,tl);
dTdx_Axy(:,:,2:end,:)   = 0.5*(dTdx(:,:,1:end-1,:) + dTdx(:,:,2:end,:));

%% dCdy on each interface
% Average value in between T-grids, in y-direction, on Axz-interface:
Ty_inbetween                = NaN(xl,yl+1,zl,tl); % One extra, to calculate difference.
Ty_inbetween(:,2:end,:,:)   = T_Axz;
%%%% Ty_inbetween(:,2:end-1,:,:) = 0.5*(tracer(:,1:end-1,:,:) + tracer(:,2:end,:,:));

% Not to waste gradients on a surface if we do not have to, we approximate
% the gradients over the first and last one using: 
Ty_inbetween(:,1,:,:)       = tracer(:,1,:,:);
%%%% Ty_inbetween(:,end,:,:)     = tracer(:,end,:,:);
%T_Axz                       = Ty_inbetween(:,2:end,:,:);
% Resulting Tracer Gradient on T-grid
dTdy = (Ty_inbetween(:,2:end,:,:) - Ty_inbetween(:,1:end-1,:,:))./repmat(dy,[1,1,1,tl]); clear Ty_inbetween

% y-Gradient on Ayz Interface, in x-direction:
dTdy_Ayz                    = 0.5*(dTdy + dTdy([2:end,1],:,:,:)); 

% y-Gradient on Axz Interface, in y-direction:
% Use original tracer values. Use that dy is constant.
dTdy_Axz                    = NaN(xl,yl,zl,tl);
dTdy_Axz(:,1:end-1,:,:)     = tracer(:,2:end,:,:) - tracer(:,1:end-1,:,:);
dTdy_Axz                    = dTdy_Axz./repmat(dy,[1,1,1,tl]);

% y-Gradient on Axy Interface, in z-direction:
dTdy_Axy                    = NaN(xl,yl,zl,tl);
dTdy_Axy(:,:,2:end,:)       = 0.5*(dTdy(:,:,1:end-1,:) + dTdy(:,:,2:end,:));

%% dCdz on each interface
% Average value in between T-grids, in z-direction, on Axy-interface:
Tz_inbetween                = NaN(xl,yl,zl+1,tl); % One extra, to calculate difference.
Tz_inbetween(:,:,2:end,:)   = T_Axy;
%%%% Tz_inbetween(:,:,2:end-1,:) = 0.5*(tracer(:,:,1:end-1,:) + tracer(:,:,2:end,:));
% Not to waste gradients on a surface if we do not have to, we approximate
% the gradients over the first and last one using: 
Tz_inbetween(:,:,1,:)       = tracer(:,:,1,:);
%%%% Tz_inbetween(:,:,end,:)     = tracer(:,:,end,:);
%T_Axy                       = Tz_inbetween(:,:,2:end,:); % Tracer value on top-Axy surfaces, 

% Resulting Tracer Gradient on T-grid.
dTdz = (Tz_inbetween(:,:,1:end-1,:) - Tz_inbetween(:,:,2:end,:))./repmat(dz,[1,1,1,tl]); clear Tz_inbetween

% z-Gradient on Ayz Interface, in x-direction:
dTdz_Ayz                    = 0.5*(dTdz + dTdz([2:end,1],:,:,:)); 

% z-Gradient on Axz Interface, in y-direction:
dTdz_Axz                    = NaN(xl,yl,zl,tl);
dTdz_Axz(:,1:end-1,:)     = 0.5*(dTdz(:,1:end-1,:) + dTdz(:,2:end,:));

% z-Gradient on Axy Interface, in z-direction:
% Use original Tracer values.
%    NOTE, we want dz=surface-bottom. But diff(z) = z2-z1 = bottom-surface.
%    We thus need -dz. BUT! our z is directed positive-up. Thus we SHOULD
%    have done diff(-z) while we did diff(z). This gives another minus sign.
%    Thus we use simply dz :).
%    Note we must use dz between tracer locations and not between edges.
dz_C1(1,1,:)                = diff(zt); % Depth difference between two t-grid interfaces. diff(zt) = z2-z1 = Bottom-Surface             
dz_C2                       = repmat(dz_C1,[xl,yl,1]);        % Make lon x lat x depth = 100x100x19
dTdz_Axy                    = NaN(xl,yl,zl,tl);
dTdz_Axy(:,:,2:end,:)       = (tracer(:,:,1:end-1,:) - tracer(:,:,2:end,:))./repmat(dz_C2,[1,1,1,tl]);
clear dz_C2 dz2

%% REMOVE SPIKES
% This part removes the spikes from the gradients.
% This is done by removing a certain percentile of the data that has values
% way beyond the average values.
% Ofcourse there are some rather arbritrary choices made here.
% if Spike==1
%    [dTdx,~,~,~,~,~,~]       = remove_999_percentale_of_tracer_gradients(dTdx);
%    [dTdy,~,~,~,~,~,~]       = remove_999_percentale_of_tracer_gradients(dTdy);
%    [dTdz,~,~,~,~,~,~]       = remove_999_percentale_of_tracer_gradients(dTdz);
%    [dTdx_Ayz,~,~,~,~,~,~]   = remove_999_percentale_of_tracer_gradients(dTdx_Ayz);
%    [dTdx_Axz,~,~,~,~,~,~]   = remove_999_percentale_of_tracer_gradients(dTdx_Axz);
%    [dTdx_Axy,~,~,~,~,~,~]   = remove_999_percentale_of_tracer_gradients(dTdx_Axy);
%    [dTdy_Ayz,~,~,~,~,~,~]   = remove_999_percentale_of_tracer_gradients(dTdy_Ayz);
%    [dTdy_Axz,~,~,~,~,~,~]   = remove_999_percentale_of_tracer_gradients(dTdy_Axz);
%    [dTdy_Axy,~,~,~,~,~,~]   = remove_999_percentale_of_tracer_gradients(dTdy_Axy);
%    [dTdz_Ayz,~,~,~,~,~,~]   = remove_999_percentale_of_tracer_gradients(dTdz_Ayz);
%    [dTdz_Axz,~,~,~,~,~,~]   = remove_999_percentale_of_tracer_gradients(dTdz_Axz);
%    [dTdz_Axy,~,~,~,~,~,~]   = remove_999_percentale_of_tracer_gradients(dTdz_Axy);
% end