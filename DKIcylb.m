function S=DKIcylb(par,b)

S0=par(1);
theta=par(2);
phi=par(3);
Dl=par(4);
Dt=par(5);
Wl=par(6);
Wt=par(7);
Wbar=par(8);


[nb]=size(b,3);
u = [sin(theta)*cos(phi),sin(theta)*sin(phi),cos(theta)]';

Trb=zeros(1,nb);
ubu=zeros(1,nb);
ubbu=zeros(1,nb);
Trbb=zeros(1,nb);
for i=1:nb
    Trb(i)=trace(b(:,:,i));
    ubu(i)=u'*b(:,:,i)*u;
    Trbb(i)=sum(sum(b(:,:,i).^2));
    ubbu(i)=u'*b(:,:,i)*b(:,:,i)*u;
end

Wbb=1/2*(10*Wt+5*Wl-15*Wbar)*ubu.^2+1/2*(5*Wbar-Wl-4*Wt)*(ubu.*Trb+2*ubbu)+Wt/3*(Trb.^2+2*Trbb);
Db=Trb*Dt+(Dl-Dt)*ubu;


Dbar=(2*Dt+Dl)/3;

S=S0*exp(-Db+ 1/6*Dbar.^2.*Wbb)';


end

