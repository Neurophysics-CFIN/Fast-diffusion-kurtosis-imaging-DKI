function [S,J] = DTI(pars,bmats)
% DTI model with parameters: S0, three Euler angles (XZX), three diffusion
% tensor eigenvalues.
S0 = pars(1);
phi = pars(2);
theta = pars(3);
psi = pars(4);
D1 = abs(pars(5));
D2 = abs(pars(6));
D3 = abs(pars(7));

R1 = [1 0 0; 0 cos(phi) sin(phi); 0 -sin(phi) cos(phi)];
R2 = [cos(theta) -sin(theta) 0;sin(theta) cos(theta) 0; 0 0 1];
R3 = [1 0 0; 0 cos(psi) sin(psi); 0 -sin(psi) cos(psi)];
R = R1*R2*R3;
RDR = R'*diag([D1 D2 D3])*R;
RDR = RDR(:);

bmats = reshape(bmats,9,[]);
nb = size(bmats,2);

S = S0*exp(-dot(bmats,RDR(:,ones(nb,1))))';

% calculate jacobian if required
if nargout==2
    J = zeros(length(S),length(pars));
    
    J(:,1) = S/S0;    
    
    dR1 = [0 0 0; 0 -sin(phi) cos(phi); 0 -cos(phi) -sin(phi)]*R2*R3;
    dR2 = R1*[-sin(theta) -cos(theta) 0; cos(theta) -sin(theta) 0; 0 0 0]*R3;
    dR3 = R1*R2*[0 0 0; 0 -sin(psi) cos(psi); 0 -cos(psi) -sin(psi)];
    dDa1 = dR1'*diag([D1 D2 D3])*R;
    dDa2 = dR2'*diag([D1 D2 D3])*R;
    dDa3 = dR3'*diag([D1 D2 D3])*R;
    dDa1 = dDa1 + dDa1';
    dDa2 = dDa2 + dDa2';
    dDa3 = dDa3 + dDa3';
    dDa1 = dDa1(:);
    dDa2 = dDa2(:);
    dDa3 = dDa3(:);
    J(:,2) = -dot(bmats,dDa1(:,ones(nb,1)))'.*S;
    J(:,3) = -dot(bmats,dDa2(:,ones(nb,1)))'.*S;
    J(:,4) = -dot(bmats,dDa3(:,ones(nb,1)))'.*S;
    
    Rk1 = R(1,:)'*R(1,:);
    Rk2 = R(2,:)'*R(2,:);
    Rk3 = R(3,:)'*R(3,:);
    Rk1 = sign(pars(5))*Rk1(:);
    Rk2 = sign(pars(6))*Rk2(:);
    Rk3 = sign(pars(7))*Rk3(:);
    J(:,5) = -dot(bmats,Rk1(:,ones(nb,1)))'.*S;
    J(:,6) = -dot(bmats,Rk2(:,ones(nb,1)))'.*S;
    J(:,7) = -dot(bmats,Rk3(:,ones(nb,1)))'.*S;
end

