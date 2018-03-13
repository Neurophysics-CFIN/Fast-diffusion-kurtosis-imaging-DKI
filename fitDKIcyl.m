function [pars,sse] = fitDKIcyl(image,B,mask,eps,imax)
% Axial DKI fitting pipeline implemented by Jonas Olesen.
%
%
% INPUT:
% image - Array of image data to be fitted to.  The first dimensions must
% discriminate between pixels, while the last dimension should correspond
% to different B-matrices. Thus image data could for instance be structured
% as [X,Y,Z,N] or [X,Y,N].
% 
% B - B-matrices having dimension 3x3xN.
%
% mask - logical array specifying which pixels to include -- if image is
% [X,Y,Z,N] then mask is [X,Y,Z]. Can optionally be left unspecified in
% which case every pixel is included.
%
% eps - Desired precission. Convergence according to step size. Smaller 
% values gives better precission but longer convergence time. This 
% parameter is optional and defaults to 1e-4.
%
% imax - Maximum number of iterations of fitting algorithm. This 
% parameter is optional and defaults to 1e4.
%
%
%
% OUTPUT:
% pars - Array of M fitting parameters arranged with the parameters as last
% index and pixels as the dist indices. If image is [X,Y,Z,N], pars is
% [X,Y,Z,M]. The order of the fitting parameters is given below.
% 
% sse - Array of sum of squared errors. If image is [X,Y,Z,N], sse is
% [X,Y,Z].
%
%
% Fitting parameters: S0, theta, phi, Dpar, Dper, Wpar, Wper, Wbar.

if nargin<4 || numel(eps)==0
    eps = 1e-4;
end
if nargin<5 || numel(imax)==0
    imax = 1e4;
end


%% adjust image dimensions and assert
dimsOld = size(image);
[image,mask] = imageAssert(image,mask);

dims = size(image);
N = prod(dims(1:end-1));
image = reshape(image,N,[]);
mask = mask(:);


%% make fit
Nfit = sum(mask);
dataFit = image(mask,:);

% initial DTI fit
DTIpars = fitDTI(dataFit,B,[],eps,imax)';


% for parallel workers (to not make DTIpars a broadcast variable)
S0s = DTIpars(1,:);
phis = DTIpars(2,:);
thetas = DTIpars(3,:);
psis = DTIpars(4,:);
Darray = DTIpars(5:7,:);

% DKI fit
parsFit = zeros(8,Nfit);
sseFit = zeros(Nfit,1);
fitfun = @(pars) DKIcyl(pars,B);
lambda = 1e-2;
upFactor = 5;
downFactor = 3;
parfor i = 1:Nfit
    par0 = ones(8,1);
    par0(1) = S0s(i);
    phi = phis(i);
    theta = thetas(i);
    psi = psis(i);
    R1 = [1 0 0; 0 cos(phi) sin(phi); 0 -sin(phi) cos(phi)];
    R2 = [cos(theta) -sin(theta) 0;sin(theta) cos(theta) 0; 0 0 1];
    R3 = [1 0 0; 0 cos(psi) sin(psi); 0 -sin(psi) cos(psi)];
    R = R1*R2*R3;
    D = Darray(:,i);
    [D,I] = sort(D,'descend');
    R = R(I,:);
    v = R(1,:)';
    par0(2) = acos(v(3));
    par0(3) = atan2(v(2),v(1));
%     par0(4) = D(1);
%     par0(5) = mean(D(2:3));
    
    [parsFit(:,i),sseFit(i)] = levmar(fitfun,par0,dataFit(i,:)',lambda,upFactor,downFactor,0,imax,eps);
end

pars = zeros(8,N);
pars(:,mask) = parsFit;
sse = zeros(N,1);
sse(mask) = sseFit;


%% adjust output to match input dimensions
if ~(numel(dimsOld)==2 && dimsOld(2)==1)
    pars = reshape(permute(pars,[2 1]),[dimsOld(1:end-1) 8]);
    sse = reshape(sse,[dimsOld(1:end-1) 1]);
end

