%Authors: Bensingh Dhas and Md Masiur Rahaman
%Description: Constructs the boundary matrix for Lagrange multiplier.
%mesh: is the mesh object containing nodes and element connectivity.
%bc: is the boundary condition object containing description of lagrange
%multipliers
%The boundary matrix is of the same size as the global stifness matrix.
%Created on: 15Aug, 2017

function [G] = lagrange_boundaryMatrix(mesh,bc)

    dofpn=2;
    dataGauss1D=gaussData1D();
    noDof=dofpn*mesh.noNd+size(bc.lagrangeuDof,1)+size(bc.lagrangevDof,1);
    G = spalloc(noDof,noDof,mesh.noNd);

    for itEg=1:size(bc.dispEdge,1)

        consNodes=bc.dispEdge(itEg,1:2);
        edgeDofu=consNodes(:);
        edgeDofv=mesh.noNd+consNodes(:);

        % Construct linear form due to applied displacement
        edgeNodes=mesh.nodes(consNodes,:);

        % Construct Lagrange multiplier boundary matrix for an edge
       
        ge = edgeLagrangeMutiplier(edgeNodes,dataGauss1D.Pt,dataGauss1D.Wt);

        % Assemble the lagrange multiplier matrix
        currEdge=bc.dispEdge(itEg,1:2);
        if bc.dispEdge(itEg,3)==1

            edgeLagDof(1)=bc.lagrangeuDof(find(bc.lagrangeuDof(:,1)==currEdge(1)),2);
            edgeLagDof(2)=bc.lagrangeuDof(find(bc.lagrangeuDof(:,1)==currEdge(2)),2);
            G(edgeDofu,edgeLagDof)= G(edgeDofu,edgeLagDof)+ge;
        end
        if bc.dispEdge(itEg,4)==1

            edgeLagDof(1)=bc.lagrangevDof(find(bc.lagrangevDof(:,1)==currEdge(1)),2);
            edgeLagDof(2)=bc.lagrangevDof(find(bc.lagrangevDof(:,1)==currEdge(2)),2);
            G(edgeDofv,edgeLagDof)= G(edgeDofv,edgeLagDof)+ge;
        end

    end
    G=G+G';
end
