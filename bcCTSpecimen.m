%Authors: Bensingh Dhas and Md Masiur Rahaman
%Description: Description of the boundary conditions for the ct-specimen.
%No boundary condition applied on equivalent strain rate.
%Created on: 15Aug, 2017
classdef bcCTSpecimen< handle
    properties
        lagrangeuDof
        lagrangevDof
        dispEdge
        noDof
        noLagDof 
        trackNode
        bcDof
    end
    methods
        function obj= bcCTSpecimen(me)
            dofpn=2;
            temp=find(me.nodes(:,1)<0.001);
            leftNodes=me.nodes(temp,:);
            temp1=intersect(find(leftNodes(:,2)>0.016), find(leftNodes(:,2)<0.019));
            temp2=intersect(find(leftNodes(:,2)<-0.016), find(leftNodes(:,2)>-0.019));
            topNodes=temp(temp1,:);
            bottomNodes=temp(temp2,:);
            obj.dispEdge=[topNodes(:)' 1 1 0 1; bottomNodes(:)' 1 1 0 -1];
            
            % lagrange multiplier table
            laguDof=[topNodes(:);bottomNodes(:)];
            obj.lagrangeuDof=[laguDof 2*me.noNd+(1:size(laguDof,1))'];
            lagvDof=[topNodes(:);bottomNodes(:)];
            obj.lagrangevDof=[lagvDof 2*me.noNd+size(laguDof,1)+(1:size(lagvDof,1))'];
            
            noLaguDof=numel(laguDof);
            noLagvDof=numel(lagvDof);
            obj.noLagDof=noLaguDof+noLagvDof;
            obj.noDof=dofpn*me.noNd+obj.noLagDof;
            obj.trackNode=topNodes(2);
            obj.bcDof=[topNodes(:); bottomNodes(:); me.noNd+topNodes(:);...
                           me.noNd+bottomNodes(:)];
            
            %Display Status
            msg=strcat('No DOF including lagrange multiplier: ', num2str(obj.noDof));
            disp(msg);
        end
    end
end