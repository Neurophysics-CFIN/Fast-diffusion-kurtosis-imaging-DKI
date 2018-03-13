function [S,J] = DKIcyl(pars,bmats)
% Axial DKI model with parameters: S0, theta, phi, Dpar, Dper, Wpar, Wper,
% Wbar.
S0 = pars(1);
theta = pars(2);
phi = pars(3);
Dpar = pars(4);
Dper = pars(5);
Wpar = pars(6);
Wper = pars(7);
Wbar = pars(8);

u = [sin(theta)*cos(phi) sin(theta)*sin(phi) cos(theta)]';

nb = size(bmats,3);
I = logical(eye(3));
Trb = sum(reshape(bmats(I(:,:,ones(nb,1))),[3 nb]))';
ubu = sum(sum(bsxfun(@times,bmats,u*u')),2);
ubu = ubu(:);
Trbb = sum(sum(bmats.^2),2);
Trbb = Trbb(:);
ubbu = sum(sum(bsxfun(@times,bmats,u(:,ones(3,1)))).^2,2);
ubbu = ubbu(:);

Wbb = 1/2*(10*Wper+5*Wpar-15*Wbar)*ubu.^2+1/2*(5*Wbar-Wpar-4*Wper)*(ubu.*Trb+2*ubbu)+Wper/3*(Trb.^2+2*Trbb);
Db = Trb*Dper+(Dpar-Dper)*ubu;
Dbar = (2*Dper+Dpar)/3;

S = S0*exp(-Db+1/6*Dbar^2.*Wbb);


if nargout==2
    J = zeros(nb,length(pars));
    
    J(:,1) = S/S0;
    
    W1 = (10*Wper+5*Wpar-15*Wbar);
    W2 = (5*Wbar-Wpar-4*Wper);
    DuDtheta = [cos(theta)*cos(phi) cos(theta)*sin(phi) -sin(theta)]';
    DubuDtheta = 2*sum(sum(bsxfun(@times,bmats,u*DuDtheta')),2);
    DubuDtheta = DubuDtheta(:);
    DubbuDtheta = 2*sum(sum(bsxfun(@times,bmats,u(:,ones(3,1)))).*sum(bsxfun(@times,bmats,DuDtheta(:,ones(3,1)))),2);
    DubbuDtheta = DubbuDtheta(:);
    DuDphi = [-sin(theta)*sin(phi) sin(theta)*cos(phi) 0]';
    DubuDphi = 2*sum(sum(bsxfun(@times,bmats,u*DuDphi')),2);
    DubuDphi = DubuDphi(:);
    DubbuDphi = 2*sum(sum(bsxfun(@times,bmats,u(:,ones(3,1)))).*sum(bsxfun(@times,bmats,DuDphi(:,ones(3,1)))),2);
    DubbuDphi = DubbuDphi(:);
    J(:,2) = (-(Dpar-Dper)*DubuDtheta + 1/6*Dbar^2*(W1*ubu.*DubuDtheta+1/2*W2*(Trb.*DubuDtheta+2*DubbuDtheta))).*S;
    J(:,3) = (-(Dpar-Dper)*DubuDphi + 1/6*Dbar^2*(W1*ubu.*DubuDphi+1/2*W2*(Trb.*DubuDphi+2*DubbuDphi))).*S;
    
    J(:,4) = (-ubu+Dbar/9*Wbb).*S;
    J(:,5) = (-Trb+ubu+2*Dbar/9*Wbb).*S;
    
    J(:,6) = (Dbar^2/6*(5/2*ubu.^2-1/2*(ubu.*Trb+2*ubbu))).*S;
    J(:,7) = (Dbar^2/6*(5*ubu.^2-2*(ubu.*Trb+2*ubbu)+1/3*(Trb.^2+2*Trbb))).*S;
    J(:,8) = (Dbar^2/6*(-15/2*ubu.^2+5/2*(ubu.*Trb+2*ubbu))).*S;
end

