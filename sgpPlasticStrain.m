%Authors: Md Masiur Rahaman and Bensingh Dhas 
%Description: Solves the gurtin mode of visco-plasticity.
%Load control with micro inertia.
%Microforce balance as the outer loop.
%Created on: 30 th Aug, 2017

function gpEpsilonp= sgpPlasticStrain(elNodes,elDisp0,elGammaDot0,elGammaDot,gpEpsilonp0,dt,dataGauss,mat)

D = mat.D; nu = mat.nu;

for itGp=1:numel(dataGauss.Wt)
    x=dataGauss.Pt(itGp,:);
    [sp,~]=Q4ShapeFn(x(1),x(2));
    B= Q4StrainDisplacement(x(1),x(2),elNodes);
    elastic_strain=B*elDisp0(:)-gpEpsilonp0(itGp,:)';
    stress=D*elastic_strain;
    
    sig_3D = [stress(1)  stress(3) 0;  stress(3) stress(2) 0 ; 0  0  nu*(stress(1)+stress(2))];
    sig_3D_dev = sig_3D - (trace(sig_3D)/3)*eye(3,3);
    tau = sqrt(sum(sum(sig_3D_dev.^2)));
    stress_dev_vec=[sig_3D_dev(1,1); sig_3D_dev(2,2); sig_3D_dev(1,2)];
    
    if tau>1e-30
        unitDevStress= stress_dev_vec/tau;
    else
        unitDevStress = zeros(numel(stress_dev_vec),1);
    end
    
    %updating gauss point epsilon p:
    
    gpGammaDot0=sp*elGammaDot0(:);
    gpGammaDot=sp*elGammaDot(:);
    gpEpsilonp(itGp,:)=gpEpsilonp0(itGp,:)+(dt/2)*(gpGammaDot0 + gpGammaDot)*unitDevStress';
    
end

end


