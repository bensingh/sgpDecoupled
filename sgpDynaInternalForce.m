%Authors: Bensingh Dhas and Md Masiur Rahaman
%Description: Computes the nodal force Q4 elements, for dynamic case. Uses
%plastic strain to determine internal forces.
%Dependencies: None
%Created on: 10th Dec, 2016

function fe=sgpDynaInternalForce(elNodes,elDisp,elDisPdDot,elEpsilonP,dataGauss,mat)
D=mat.D; rho = mat.rho;
fe=zeros(8,1);
for itGp=1:numel(dataGauss.Wt)
    x=dataGauss.Pt(itGp,:);
    B= Q4StrainDisplacement(x(1),x(2),elNodes);
    [sp,dsp]=Q4ShapeFn(x(1),x(2));
    JxW=det(dsp*elNodes)*dataGauss.Wt(itGp);
    stress=D*(B*elDisp(:)-elEpsilonP(itGp,:)');
    fe(1:8,1)=fe(1:8,1)+JxW*B'*stress;
    
    fe(1:4,1)=fe(1:4,1)+JxW*rho*(sp'*sp)*elDisPdDot(:,1);
    fe(5:8,1)=fe(5:8,1)+JxW*rho*(sp'*sp)*elDisPdDot(:,2);
end
end