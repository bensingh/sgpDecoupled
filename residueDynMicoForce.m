%Authors: Bensingh Dhas
%Description: Constructs the residue vector associated with
%microforce equation of gurtin includes micro inertia.
%Created on: 17th Jan, 2017.

function fMf=residueDynMicoForce(mesh,gammaP,gammaPdot,gammaPdDot,gp_stress,dataGauss,mat)
  fMf=zeros(mesh.noNd,1);
  nu=mat.nu;
  s0=mat.S_0;l1=mat.l1;
  h0=mat.H; d_0=mat.d_0;
  m1=mat.m1;
  lm = mat.lm ; rho = mat.rho;
  for it_el=1:mesh.noEl
    nodeIndex=mesh.elements(it_el,1:4);
    el_nodes=mesh.nodes(nodeIndex,:);
    el_gammaP=gammaP(nodeIndex);
    el_gammaPdot=gammaPdot(nodeIndex,1);
    el_gammaPdDot=gammaPdDot(nodeIndex,1);
    el_stress=squeeze(gp_stress(it_el,:,:));
    for itGp=1:numel(dataGauss.Wt)
      x=dataGauss.Pt(itGp,:);
      [sp,dsp]=Q4ShapeFn(x(1),x(2));
      jacobian=dsp*el_nodes;
      JxW=det(jacobian)*dataGauss.Wt(itGp);
      dsp=jacobian\dsp;

      gpGammaP=sp*el_gammaP(:);
      gpGammaPdot=sp*el_gammaPdot(:);
      gpGammaPdot(gpGammaPdot<1e-10)=1e-10;
      gpGammaPdDot=sp*el_gammaPdDot(:);
    
    stress = el_stress(itGp,:)';
    sig_3D = [stress(1)  stress(3) 0;  stress(3) stress(2) 0 ; 0  0  nu*(stress(1)+stress(2))];
    sig_3D_dev = sig_3D - (trace(sig_3D)/3)*eye(3,3);
    tau = sqrt(sum(sum(sig_3D_dev.^2)));
    
    f4= (1+(h0/s0)*gpGammaP)*((gpGammaPdot/d_0)^m1);
    f5_dyn =(rho/s0)*lm^2*gpGammaPdDot;
    
    fMf(nodeIndex,1)=fMf(nodeIndex,1)+JxW*l1^2*(dsp'*dsp*el_gammaP(:));
    fMf(nodeIndex,1)=fMf(nodeIndex,1)+JxW*sp'*(f4+f5_dyn-(tau/s0));

    end
  end
end