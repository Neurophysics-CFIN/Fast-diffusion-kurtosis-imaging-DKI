function metrics = DKIcylMetrics(DKIcylPars)
% Calculates DKI metrics from axially symmetric DKI fitting parameters.
%
% INPUT:
% DKIcylPars - fitting parameters from axially symmetric DKI fit. The first
% three parameters are not used but the following parameters must be in the
% order Dpar, Dper, Wpar, Wper, and Wbar. The parameters must be contained
% in the last index while the first indices corresponds to pixels -- for
% instance [X,Y,Z,N] where N=8.
%
% OUTPUT
% metrics - Contains the DKI metrics Dpar, Dper, Dbar, Wpar, Wper, Wbar,
% Kpar, Kper, and Kbar in that order. These are contained in the last index
% while the first indices corresponds to pixels and are in the same format
% as the input.

dimsOld = size(DKIcylPars);
DKIcylPars = imageAssert(DKIcylPars,[]);

Dpar = DKIcylPars(:,:,:,4);
Dper = DKIcylPars(:,:,:,5);
Wpar = DKIcylPars(:,:,:,6);
Wper = DKIcylPars(:,:,:,7);
Wbar = DKIcylPars(:,:,:,8);

Dbar = (Dpar+2*Dper)/3;
Kpar = Wpar.*Dbar.^2./Dpar.^2;
Kper = Wper.*Dbar.^2./Dpar.^2;
Kbar = 2/315*(Dpar.^2.*(75*Wbar+20*Wpar-32*Wper)+4*Dpar.*Dper.*(30*Wbar-Wpar-8*Wper)+8*Dper.^2.*(15*Wbar-2*Wpar+8*Wper))./(2*Dbar.^2);

dims = size(DKIcylPars);
metrics = zeros([dims(1:3), 9]);
metrics(:,:,:,1) = Dpar;
metrics(:,:,:,2) = Dper;
metrics(:,:,:,3) = Dbar;
metrics(:,:,:,4) = Wpar;
metrics(:,:,:,5) = Wper;
metrics(:,:,:,6) = Wbar;
metrics(:,:,:,7) = Kpar;
metrics(:,:,:,8) = Kper;
metrics(:,:,:,9) = Kbar;

metrics = reshape(metrics,[dimsOld(1:end-1) 9]);