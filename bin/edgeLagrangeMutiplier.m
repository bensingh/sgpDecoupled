%Authors: Bensingh Dhas and Md Masiur Rahaman
%Description: Computes the lagrange multiplier matrix for a Q4 edge in
%dimensional elasticity problems.
%Input: Element nodes and the edges for which u and v directions has to be
%fixed. The edges has the numbering convention of 1 to 4 counter clockwise
%with bottom edge numbered 1. The assembly is in the usual convention {u,v}
%Dependencies: None
%Created on: 07th Oct, 2016

function biliEdge=edgeLagrangeMutiplier(edgeNodes,gauss_pt,gauss_wt)

    biliEdge=zeros(2,2);
    edgeLength=norm(edgeNodes(1,:)-edgeNodes(2,:));
    % Construct Lagrange multiplier matrix for an edge
    for it_gp=1:length(gauss_pt)
        xi=gauss_pt(it_gp);
        [sp,~]=P2ShapeFn(xi);
        biliEdge(1:2,1:2)=biliEdge(1:2,1:2)+gauss_wt(it_gp)*edgeLength*sp(:)*sp(:)';
    end
     
end
