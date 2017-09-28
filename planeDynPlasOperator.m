%Authors: Bensingh Dhas
%Description: Constructs the tangent matrix associated with
%microforce equation of gurtin includes micro inertia.
%Created on: 21st Nov, 2016
function kMf=planeDynPlasOperator(mesh,gammaP,gammaPdot,dt,dataGauss,mat)

kMf=sparse(mesh.noNd,mesh.noNd,mesh.noNd);
s0=mat.S_0;l1=mat.l1;
h0=mat.H; d_0=mat.d_0;
m1=mat.m1;
lm = mat.lm ; rho = mat.rho;

for it_el=1:mesh.noEl
  nodeIndex=mesh.elements(it_el,1:4);
  el_nodes=mesh.nodes(nodeIndex,:);
  el_gammaPdot=gammaPdot(nodeIndex,1);
  el_gammaP=gammaP(nodeIndex,1);
  ke=zeros(4,4);
  for itGp=1:numel(dataGauss.Wt)
    x=dataGauss.Pt(itGp,:);
    [sp,dsp]=Q4ShapeFn(x(1),x(2));
    jacobian=dsp*el_nodes;
    JxW=det(jacobian)*dataGauss.Wt(itGp);
    dsp=jacobian\dsp;
    
    gpGammaPdot=sp*el_gammaPdot(:);
    gpGammaP=sp*el_gammaP(:);
    gpGammaPdot(gpGammaPdot<1e-10)=1e-10;
    
    c4= (h0/s0)*((gpGammaPdot/d_0)^m1);
    c5= (2/dt)*(m1/d_0^m1)*(1+(h0/s0)*gpGammaP)*gpGammaPdot^(m1-1);
    c5_dyn = (rho/s0)*lm^2*(4/dt^2);
    
    ke=ke+JxW*l1^2*(dsp'*dsp);
    ke=ke+JxW*(c4+c5+c5_dyn)*(sp'*sp);
    
  end
  kMf(nodeIndex,nodeIndex)=kMf(nodeIndex,nodeIndex)+ke;
end
end
