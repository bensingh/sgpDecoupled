%Created by Bensingh Dhas
%Description: Computes the stress at the gauss points
%of the element.
%Dependencies: None
function stress = sgpElementStress(elNodes,elDisp,elEpsilonP,dataGauss,mat)

D = mat.D;
elDisp=elDisp(:);

for it_gp=1:numel(dataGauss.Wt)
    
    x=dataGauss.Pt(it_gp,:);
    B= Q4StrainDisplacement(x(1),x(2),elNodes);
    strain=(B*elDisp)-elEpsilonP(it_gp,:)';
    stress(it_gp,:)=(D*strain)';
    
end

end
