function [T_Ayz,T_Axz,T_Axy] = tracer_on_interface_xyzt(tracer)
% This script calculAtes tracer values on the interfaces surrounding the
% T-grid. It is made for [x,y,z,t]-coordinates.

%% Pre calculations and explanations:
% Length of vectors.
[xl,yl,zl,tl]   = size(tracer);

%% dCdx on each interface
% Policy:  There is a value on each interface "east" of the T-grid point.
% Average value of tracer, on interface between two tracer-values.
% Interface is to the East, in x-direction, on Ayz-interface:
T_Ayz            = 0.5*(tracer + tracer([2:end,1],:,:,:));

%% dCdy on each interface
% Policy: There is a value on each interface "north" of the T-grid point.
% Average value in between T-grids, in y-direction, on Axz-interface:
T_Axz                = NaN(xl,yl,zl,tl);
T_Axz(:,1:end-1,:,:) = 0.5*(tracer(:,1:end-1,:,:) + tracer(:,2:end,:,:));
% Not to waste gradients on a surface if we do not have to, we approximate
% the gradients over the first and last one using: 
T_Axz(:,end,:,:)     = tracer(:,end,:,:);

%% dCdz on each interface
% Policy: There is a value on each interface "bottom" of the T-grid point.
% Average value in between T-grids, in z-direction, on Axy-interface:
T_Axy                = NaN(xl,yl,zl,tl);
T_Axy(:,:,1:end-1,:) = 0.5*(tracer(:,:,1:end-1,:) + tracer(:,:,2:end,:));
% Not to waste gradients on a surface if we do not have to, we approximate
% the gradients over the first and last one using: 
T_Axy(:,:,end,:)     = tracer(:,:,end,:);