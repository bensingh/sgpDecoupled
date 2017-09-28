%Authors: Bensingh Dhas and Md Masiur Rahaman
%Description: This routine calculates the generalised force for associated with lagrange
%multiplier alone.
%Created on: 21Aug, 2017
function  globalForce=lagrangeForce(mesh,bc,appDisp,noDof)
    globalForce=zeros(noDof,1);
    dataGauss1D=gaussData1D();
    for it_eg=1:size(bc.dispEdge,1)
      
        consNodes=bc.dispEdge(it_eg,1:2);
        
        % Construct Lagrange multiplier boundary matrix for an edge
        edgeNodes=mesh.nodes(consNodes,:);
        edgeLength=norm(edgeNodes(1,:)-edgeNodes(2,:));
        biliEdge=zeros(2,2);
        for itGp=1:numel(dataGauss1D.Wt)
            xi=dataGauss1D.Pt(itGp);
            JxW=0.5*dataGauss1D.Wt(itGp)*edgeLength;
            [sp,~]=P2ShapeFn(xi);
            biliEdge(1:2,1:2)=biliEdge(1:2,1:2)+JxW*sp(:)*sp(:)';
        end
        
        if bc.dispEdge(it_eg,3)==1
            edgeDofu=consNodes(:);
            currEdge=bc.dispEdge(it_eg,1:2);
            edgeLagDof(1)=bc.lagrangeuDof(find(bc.lagrangeuDof(:,1)==currEdge(1)),2);
            edgeLagDof(2)=bc.lagrangeuDof(find(bc.lagrangeuDof(:,1)==currEdge(2)),2);
            unBalDisp=appDisp(edgeLagDof);
            globalForce(edgeLagDof,1)= globalForce(edgeLagDof,1)+biliEdge*unBalDisp(:);
        end
        if bc.dispEdge(it_eg,4)==1
            edgeDofv=mesh.noNd+consNodes(:);
            currEdge=bc.dispEdge(it_eg,1:2);
            edgeLagDof(1)=bc.lagrangevDof(find(bc.lagrangevDof(:,1)==currEdge(1)),2);
            edgeLagDof(2)=bc.lagrangevDof(find(bc.lagrangevDof(:,1)==currEdge(2)),2);
            unBalDisp=appDisp(edgeLagDof);
            globalForce(edgeLagDof,1)= globalForce(edgeLagDof,1)+biliEdge*unBalDisp(:);
        end

    end
end
