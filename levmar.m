function [x,s2,n] = levmar(f,x,y,lambda,upFactor,downFactor,gradeps,imax,eps)
% Levenberg-Marquardt implementation by Jonas Olesen for unbounded
% non-linear fitting based on the algortihm presented by Transtrum and
% Sethna.
% 
%
% INPUT:
% f - Objective function in the format [y,J] = f(x), where y contains the
% function values at x and J is the Jacobian.
%
% x - Initial point to start algorithm from.
% 
% y - Function values to fit to.
% 
% lambda - Parameter for weighting between Gauss-Newton and
% gradient-descent.
% 
% upFactor/downFactor - Scaling factors for running adjustment of lambda.
% 
% gradeps - Desired precission. Convergence according to gradient. Smaller 
% values gives better precission but longer convergence time. Can be set to
% zero to use "eps" exclusively.
% 
% imax - Maximum number of iterations at which the algorithm stops even if
% not converged.
%
% eps - Desired precission. Convergence according to step size. Smaller 
% values gives better precission but longer convergence time. Can be set to
% zero to use "gradeps" exclusively. This parameter is taken as optional
% and defaults to 0.
%
%
% OUTPUT:
% x - Optimal parameters found.
% 
% s2 - Sum of squared errors at optimal point.
%
% n - Number of iterations.


if nargin<9
    eps = 0;
end


%%
minDiagVals = 1e-3;
relDxTol = 1e-12;

warningstate1 = warning('off','MATLAB:nearlySingularMatrix');
warningstate2 = warning('off','MATLAB:singularMatrix');
warningstate3 = warning('off','MATLAB:rankDeficientMatrix');

maxDiagFound = minDiagVals*ones(length(x),1);

x = x(:);
[fx,J] = f(x);
r = y - fx;
s2 = r'*r;
for n = 1:imax
    Jtr = J'*r;
    if norm(2*Jtr)<gradeps
        break
    end
    
    JtJ = J'*J;
    D = diag(JtJ);
    filter = D<maxDiagFound;
    D(filter) = maxDiagFound(filter);
    maxDiagFound(~filter) = D(~filter);
    DtD = diag(D);
    
    dx = (JtJ + lambda*DtD)\Jtr;
    
    xtrial = x+dx;
    [rtrial,s2trial] = S2(f,xtrial,y);
    if s2trial<s2
        x = xtrial;
        r = rtrial;
        s2 = s2trial;
        [~,J] = f(x);
        if lambda>1e-9
            lambda = lambda/downFactor;
        end
    else
        if norm(dx./x)<relDxTol
%             fprintf('levmar: failed to reduce gradient norm below specified value\n')
            return
        end
        lambda = lambda*upFactor;
    end
    
    if norm(dx./x)/length(dx)<eps
        break
    end
end

if n==imax
    fprintf('levmar: reached maximum number of iterations\n')
end

warning(warningstate1)
warning(warningstate2)
warning(warningstate3)
end



function [r,s2] = S2(f,x,y)
fx = f(x);
r = y - fx;
s2 = r'*r;
end