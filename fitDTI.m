function [pars,sse] = fitDTI(image,B,mask,eps,imax)
% DTI fitting pipeline implemented by Jonas Olesen.
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
% imax - Maximum number of iterations of fitting algorithm.
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
% Fitting parameters (order): S0, three Euler angles (XZX), three diffusion
% tensor eigenvalues.

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

parsFit = zeros(7,Nfit);
parsFit(1,:) = dataFit(:,1);
parsFit(5:7,:) = 0.5;
sseFit = zeros(Nfit,1);
fitfun = @(pars) DTI(pars,B);
lambda = 1e-2;
upFactor = 5;
downFactor = 3;
parfor i = 1:Nfit
    [parsFit(:,i),sseFit(i)] = levmar(fitfun,parsFit(:,i),dataFit(i,:)',lambda,upFactor,downFactor,0,imax,eps);
end
parsFit(5:7,:) = abs(parsFit(5:7,:));

pars = zeros(7,N);
pars(:,mask) = parsFit;
sse = zeros(N,1);
sse(mask) = sseFit;


%% adjust output to match input dimensions
if ~(numel(dimsOld)==2 && dimsOld(2)==1)
    pars = reshape(permute(pars,[2 1]),[dimsOld(1:end-1) 7]);
    sse = reshape(sse,[dimsOld(1:end-1) 1]);
end

