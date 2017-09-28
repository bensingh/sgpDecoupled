%Authors: Bensingh Dhas and Md Masiur Rahaman
%Description: Constructs the boundary matrix for Lagrange multiplier.
%mesh: is the mesh object containing nodes and element connectivity.
%bc: is the boundary condition object containing description of lagrange
%multipliers
%The boundary matrix is of the same size as the global stifness matrix.
%Created on: 15Aug, 2017

function [G,fDisp] = boundaryMatrix(mesh,bc)

    dofpn=2;
    dataGauss1D=gaussData1D();
    noDof=dofpn*mesh.noNd+size(bc.lagrangeuDof,1)+size(bc.lagrangevDof,1);
    G = spalloc(noDof,noDof,mesh.noNd);
    fDisp = spalloc(noDof,1,mesh.noNd);
    for itEg=1:size(bc.dispEdge,1)

        consNodes=bc.dispEdge(itEg,1:2);
        edgeDofu=consNodes(:);
        edgeDofv=mesh.noNd+consNodes(:);

        % Construct linear form due to applied displacement
        edgeNodes=mesh.nodes(consNodes,:);
        edgeLength=norm(edgeNodes(1,:)-edgeNodes(2,:));

        % Construct Lagrange multiplier boundary matrix for an edge
        biliEdge=zeros(2,2);
        for itGp=1:numel(dataGauss1D.Wt)
            xi=dataGauss1D.Pt(itGp);
            JxW=dataGauss1D.Wt(itGp)*edgeLength*0.5;
            [sp,~]=P2ShapeFn(xi);
            biliEdge(1:2,1:2)=biliEdge(1:2,1:2)+JxW*sp(:)*sp(:)';
        end

        % Assemble the lagrange multiplier matrix
        currEdge=bc.dispEdge(itEg,1:2);
        if bc.dispEdge(itEg,3)==1
            udisp=bc.dispEdge(itEg,5);
            edgeLagDof(1)=bc.lagrangeuDof(find(bc.lagrangeuDof(:,1)==currEdge(1)),2);
            edgeLagDof(2)=bc.lagrangeuDof(find(bc.lagrangeuDof(:,1)==currEdge(2)),2);
            G(edgeDofu,edgeLagDof)= G(edgeDofu,edgeLagDof)+biliEdge;
            fDisp(edgeLagDof,1)=fDisp(edgeLagDof,1)+biliEdge*udisp*[1;1];
        end
        if bc.dispEdge(itEg,4)==1
            vdisp=bc.dispEdge(itEg,6);
            edgeLagDof(1)=bc.lagrangevDof(find(bc.lagrangevDof(:,1)==currEdge(1)),2);
            edgeLagDof(2)=bc.lagrangevDof(find(bc.lagrangevDof(:,1)==currEdge(2)),2);
            G(edgeDofv,edgeLagDof)= G(edgeDofv,edgeLagDof)+biliEdge;
            fDisp(edgeLagDof,1)=fDisp(edgeLagDof,1)+biliEdge*vdisp*[1;1];
        end

    end
    G=G+G';
end



