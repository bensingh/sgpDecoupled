
function integrationPath= edgeNaturalUnitNormal(pathJ)

    % Gauss Point along each edge
    bt_gp=[-sqrt(1/3),-sqrt(1/3),sqrt(1/3),-sqrt(1/3)];
    rt_gp=[sqrt(1/3),-sqrt(1/3),sqrt(1/3),sqrt(1/3)];
    tp_gp=[sqrt(1/3),sqrt(1/3),-sqrt(1/3),sqrt(1/3)];
    lt_gp=[-sqrt(1/3),sqrt(1/3),-sqrt(1/3),-sqrt(1/3)];
    edgeGaussPoint=[bt_gp;rt_gp;tp_gp;lt_gp];
    edgeUnitNormal=[1 0 -1;2 1 0;3 0 1;4 -1 0];
    integrationPath=zeros(size(pathJ,1),7);
    for i =1:size(pathJ,1)
        edNaturalNormal=edgeUnitNormal((edgeUnitNormal(:,1)==pathJ(i,2)),2:3);
        edGaussPt=edgeGaussPoint(pathJ(i,2),:);
        integrationPath(i,:)=[pathJ(i,1) edGaussPt edNaturalNormal];
    end

end

        

                
