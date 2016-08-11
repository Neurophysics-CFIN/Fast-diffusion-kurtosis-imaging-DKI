function S=DTI(par,varargin)
nargs=size(varargin,2);
if nargs==4
    g=varargin{1};
    Delta=varargin{2};
    delta=varargin{3};
    [m,n]=size(g);
    ndiff=varargin{4};
    if ~(n==3 || ((m==1 || n==1) &&  mod(numel(g),3)==0))
        disp('g wrong size')
        return;
    end
    
elseif nargs==2
    b=varargin{1};
    ndiff=varargin{2};
 
elseif nargs==1
    b=varargin{1};
    ndiff=size(b,3);
    
else
    error('Wrong number of input arguments! Check your fitscript.');
    return
end

if length(par)~=7
    disp('par wrong size');
    return;
end

S0=par(1);
phi=par(2);
theta =par(3);
psi=par(4);
D1=par(5);
D2=par(6);
D3=par(7);
D=diag([D1 D2 D3]);

gamma=2.675222000000000e-003;

R1=[1 0 0; 0 cos(phi) sin(phi); 0 -sin(phi) cos(phi)];
R2=[cos(theta) -sin(theta) 0;sin(theta) cos(theta) 0; 0 0 1];
R3=[1 0 0; 0 cos(psi) sin(psi); 0 -sin(psi) cos(psi)];

R=R1*R2*R3;

if nargs==4
    q=gamma*delta*g;
    q=reshape(q',3,[]);
    
    
    for i=1:length(q)
        y(i)=q(:,i)'*R'*D*R*q(:,i);
    end
    y=(Delta-delta/3)*y;
elseif nargs==2||nargs==1
    allb=reshape(b,9,ndiff);
    
    RDR=reshape(repmat(R'*D*R,[1,1,ndiff]),9,ndiff);
    y=dot(allb,RDR);
end
S=S0*exp(-y);
%S=S';