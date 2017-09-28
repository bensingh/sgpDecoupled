%Authors: Bensingh Dhas and Md Masiur Rahaman
% Description: Creates the System matrix for the complementary equation
%Dependencies: None
%Created on: 12th Oct, 2016

function [keOmega,feOmega]=elementComplementaryDynamic(el_nodes,el_gamma0,el_gamma_dot0,el_omegaP0,el_omegaPdot0,el_omegaPdDot0,dt,dataGauss,mat)
       
keOmega=zeros(4,4);
feOmega=zeros(4,1);
s0=mat.S_0;l1=mat.l1;
h0=mat.H; d_0=mat.d_0;
m1=mat.m1;
lm = mat.lm ; rho = mat.rho;

for it_gp=1:numel(dataGauss.Wt)

x=dataGauss.Pt(it_gp,:);
wt=dataGauss.Wt(it_gp);
[sp,dsp]=Q4ShapeFn(x(1),x(2));
jacobian=dsp*el_nodes;
det_jac=det(jacobian);
inv_jacobian=inv(jacobian);
dsp = inv_jacobian*dsp;

gp_gammaP=sp*el_gamma0(:);
gp_gammaPdot=sp*el_gamma_dot0(:);

gp_omegaP0 = sp*el_omegaP0(:);
gp_omegaPdot0 = sp*el_omegaPdot0(:);
gp_omegaPdDot0 = sp*el_omegaPdDot0(:);

gp_gammaPdot(gp_gammaPdot<1e-10)=1e-10;


c4= (h0/s0)*((gp_gammaPdot/d_0)^m1);
c5= (1/s0)*(2/dt)*(m1/d_0^m1)*(s0+h0*gp_gammaP)*gp_gammaPdot^(m1-1);
c5_dyn = (rho/s0)*lm^2*(4/dt^2);

c6_dyn = (rho/s0)*lm^2;

for i=1:4
    for j=1:4
        keOmega(i,j)=keOmega(i,j)+l1^2*(det_jac*wt*dsp(:,j)'*dsp(:,i));
        keOmega(i,j)=keOmega(i,j)+(det_jac*wt)*(c4-c5+c5_dyn)*sp(i)*sp(j);
    end
    
    feOmega(i,1)=feOmega(i,1)-(det_jac*wt)*sp(i)*c5*(gp_omegaP0+(dt/2)*gp_omegaPdot0);
    feOmega(i,1)=feOmega(i,1)+(det_jac*wt)*sp(i)*(c5_dyn*(gp_omegaP0+dt*gp_omegaPdot0)+c6_dyn*gp_omegaPdDot0);
    
end
end

end
