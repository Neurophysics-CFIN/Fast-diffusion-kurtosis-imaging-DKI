function DKIpar = FitAxDKI(img,b,mask,options,v1)
% img: image data, which is given as a RxCxN or RxCxSxN array, where R and
% C are the number of rows and colums in the images, S is the number of
% slices, and N is the number of images. b: b matrix, which is given as a
% 3x3xN array. mask: is an optional RxCxS array, which contains ones at the
% voxels to be fitted and zeros elsewhere. Alternatively, it can be given
% as []. A fit is then made to all voxels. options: is an optional argument
% containing fitting options for the lsqcurvefit function provided by
% Matlab. If no input is given standard options will be used. v1: is an
% optional argument, which is given as a RxCx3 or RxCxSx3 array depending
% on the dimensions of the img input argument. It contains the principal
% eigenvector of the diffusion tensor for each voxel. If no input is given
% a DTI fit will be carried out to determine the diffusion tensor. DKIpar:
% is the output given as a RxCx8 or RxCxSx8 array depending on the
% dimensions of the img input argument. It contains the DKI fit result for
% each voxel in the following order: non diffusion weighted signal
% amplitude, angle theta, angle phi, axial diffusivity, radial diffusivity,
% axial kurtosis, radial kurtosis, mean kurtosis. Use of this script is
% free, but please cite Hansen, Shemesh, and Jespersen, NeuroImage 2016.

% the image data is reshaped to fit the fitting procedure used below
dim = size(img);
if(length(dim)==3)
    N = dim(3);
    ni = dim(1);
    nj = dim(2);
    % setting mask to everything if unspecified
    if nargin<3 || numel(mask)==0
        mask = true(ni,nj);
    end
elseif(length(dim)==4)
    N = dim(4);
    ni = dim(1);
    nj = dim(2)*dim(3);
    img = reshape(img,[ni nj N]);
    if nargin == 5, v1 = reshape(v1,[ni nj 3]); end
    % setting mask to everything if unspecified
    if nargin<3 || numel(mask)==0
        mask = true(ni,nj);
    end
    mask = repmat(mask,[dim(3) 1]);
else
    disp('incorrect dimensions of image data!');
    return
end


% Setting unspecified input arguments.
S0 = img(:,:,1);
if nargin<4 || numel(options)==0
    options=optimset('lsqcurvefit');
    options=optimset(options,'Display','off','algorithm','levenberg-marquardt');
    options=optimset(options,'MaxIter',10000,'TolFun',1e-6,'TolX',1e-6);
    lb = []; ub = [];
    
    if nargin<5
        % DTI fit to get principal axis
        DTIfun = @(par,b) DTI(par,b);
        DTIpar = zeros(ni,nj,7);
        v1 = zeros(ni,nj,3);
        progress0 = 0;
        total = sum(mask(:));
        count = 0;
        for i = 1:ni
        for j = 1:nj
            if mask(i,j)
                y = squeeze(double(img(i,j,:)));
                par0 = zeros(7,1);
                par0(1) = S0(i,j);
                par0(4:5)=[2,0.3]';
                DTIpar(i,j,:) = lsqcurvefit(DTIfun,par0,b,y,lb,ub,options);
                
                phi = DTIpar(i,j,2);
                theta = DTIpar(i,j,3);
                psi = DTIpar(i,j,4);
                R1 = [1 0 0; 0 cos(phi) sin(phi); 0 -sin(phi) cos(phi)];
                R2 = [cos(theta) -sin(theta) 0;sin(theta) cos(theta) 0; 0 0 1];
                R3 = [1 0 0; 0 cos(psi) sin(psi); 0 -sin(psi) cos(psi)];
                R = R1*R2*R3;
        
                D = DTIpar(i,j,5:7);
                [~,I] = sort(D,'descend');
                R = R(I,:);
                v1(i,j,:) = R(1,:);
                
                count = count+1;
                progress = floor(count/total*10);
                if progress>progress0
                    progress0 = progress;
                    disp(['DTI progress: ' num2str(progress*10) '%']);
                end
            end
        end
        end
    end
end


% axially symmetric DKI fit
DKIfun = @(par,b) DKIcylb(par,b);
DKIpar = zeros(ni,nj,8);
progress0 = 0;
total = sum(sum(mask));
count = 0;
for i = 1:1:ni
for j = 1:1:nj
    if mask(i,j)
        y = squeeze(double(img(i,j,:)));
        par0 = zeros(8,1);
        par0(1) = S0(i,j);
        par0(2) = acos(v1(i,j,3));
        par0(3) = atan2(v1(i,j,2),v1(i,j,1));
        DKIpar(i,j,:) = lsqcurvefit(DKIfun,par0,b,y,lb,ub,options);
        
        count = count+1;
        progress = floor(count/total*10);
        if progress>progress0
            progress0 = progress;
            disp(['DKI progress: ' num2str(progress*10) '%']);
        end
    end
end
end


% reshape DKIpar to match input dimensions
if length(dim)==4
    DKIpar = reshape(DKIpar,[dim(1:3) 8]); % test
end
end