%Created by Bensingh Dhas
%Description: Computes the gauss points. Mostly used for plotting the stresses.
%Dependencies: None
%Created on: 6th, Oct, 2016
function meshGp = deformedGausspt(mesh,dataGauss)
  meshGp=zeros(mesh.noEl,4,2);
  for itEl=1:mesh.noEl
    el_gauss_pt=zeros(numel(dataGauss.Wt),2);
    elNodes=mesh.nodes(mesh.elements(itEl,:),:);
    for it_gp=1:numel(dataGauss.Wt)
      xi=dataGauss.Pt(it_gp,1);
      eta=dataGauss.Pt(it_gp,2);
      [sp,~]=Q4ShapeFn(xi,eta);
      el_gauss_pt(it_gp,1)=sp*elNodes(:,1);
      el_gauss_pt(it_gp,2)=sp*elNodes(:,2);
    end
    meshGp(itEl,:,:)=el_gauss_pt;
  end
end
