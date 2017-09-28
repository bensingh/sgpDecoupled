%Authors: Bensingh Dhas and Md Masiur Rahaman
%Description: Creates the plane elastic stiffness using Q4 elements.
%Dependencies: None
%Created on: 10th Dec, 2016

function ke=Q4PlaneElasticOperator(elNodes,dataGauss,mat)

    D = mat.D;
    ke=zeros(8,8);
    for itGp=1:numel(dataGauss.Wt)
        x=dataGauss.Pt(itGp,:);
        B= Q4StrainDisplacement(x(1),x(2),elNodes);
        [~,dsp]=Q4ShapeFn(x(1),x(2));
        JxW=det(dsp*elNodes)*dataGauss.Wt(itGp);
        ke(1:8,1:8)=ke(1:8,1:8)+JxW*B'*D*B;
    end

end
