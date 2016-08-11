function [Dl,Dt,MD,Wl,Wt,Wbar,Kpar,Kperp,Kbar] = DKIcylMetrics(par)
% extracts metrics from par array. It is assumed that the input par
% contains fitting parameters from DKIcylb along the last dimension.
% Output: axial diffusivity, radial diffusivity, mean diffusivity, 
% axial tensor kurtosis, radial tensor kurtosis, mean tesnor kurtosis, 
% axial kurtosis, radial kurtosis, mean kurtosis.
% Use of this script is free, but please cite Hansen, Shemesh, and
% Jespersen, NeuroImage 2016.

% Extracting data from par based on its dimensions
N = ndims(par);
if N==3
	if size(par,1)==1 && size(par,2)==1
		Dl = par(1,1,4);
		Dt = par(1,1,5);
		Wl = par(1,1,6);
		Wt = par(1,1,7);
		Wbar = par(1,1,8);
	else
		Dl = par(:,:,4);
		Dt = par(:,:,5);
		Wl = par(:,:,6);
		Wt = par(:,:,7);
		Wbar = par(:,:,8);
	end
elseif N==4
    Dl = par(:,:,:,4);
    Dt = par(:,:,:,5);
    Wl = par(:,:,:,6);
    Wt = par(:,:,:,7);
    Wbar = par(:,:,:,8);
else
    disp('par has incorrect dimensions!');
end

% calculating metrics
MD = 1/3*(Dl+2*Dt);
Dbar = (Dl+2*Dt)/3;
Kpar = Wl.*Dbar.^2./Dl.^2;
Kperp = Wt.*MD.^2./Dt.^2;
Kbar = 2/315*(Dl.^2.*(75*Wbar+20*Wl-32*Wt)+4*Dl.*Dt.*(30*Wbar-Wl-8*Wt)+8*Dt.^2.*(15*Wbar-2*Wl+8*Wt))./(2*Dbar.^2);
