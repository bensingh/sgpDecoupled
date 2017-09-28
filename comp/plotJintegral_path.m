% Pls use this routine with caution. Highly specialised for a
% crack_ct_4.inp
%
function pathElConn =plotJintegral_path(inpfname,dataGauss)
    % inpfname='solver_crack_ct_4_path1.inp';
    [nodes,connect,unit_force_vector,dead_dof_u,dead_dof_gp,pathJ,domainJ,track_node]=readNodesConnectAbaqus(inpfname);
    % hold on
    noPathEl=size(pathJ,1);
    pathElConn=connect(pathJ(:,1),:);
    connectGraph=zeros(noPathEl,noPathEl);
    for i=1:size(pathElConn,1)
        curEl=pathElConn(i,:);
        for j=1:size(pathElConn,1)
            isec=intersect(curEl,pathElConn(j,:));
            if ~isempty(isec)
                connectGraph(i,j)=1;
            end
        end
    end
    rowSum=sum(connectGraph)';
    start=find(rowSum==2);
    newConnect=start(1);
    for i=1:noPathEl
        nextEl=find(connectGraph(:,newConnect(end)));
        nextEl=setdiff(nextEl,newConnect);
        newConnect=[newConnect;nextEl];
    end
    newPathJ=pathJ(newConnect,1);
    pathElConn=connect(newPathJ(:,1),:);
    
    for i=1:45

        currEl=pathElConn(i,:);
        nextEl=pathElConn(i+1,:);
        commEdge=intersect(currEl,nextEl);
        natEdgeCur=[find(currEl==commEdge(1)),find(currEl==commEdge(2))];
        natEdgeNoCurr=floor(sum(natEdgeCur)/2); % crude logic think of it
        natEdgeCur=unique(natEdgeCur); % get natural edge in ascending order
        natEdgeNext=[find(nextEl==commEdge(1)),find(nextEl==commEdge(2))];
        natEdgeNoNext=floor(sum(natEdgeNext)/2); % crude logic think of itE
        natEdgeNext=unique(natEdgeNext); % get natural edge in ascending order
        if isequal(natEdgeCur,natEdgeNext)
            % apply two end transposition
            temp=pathElConn(i+1,:);
            temp=[temp(end) temp(1:end-1)]; % first trans position
            temp=[temp(end) temp(1:end-1)]; % second trans position
            pathElConn(i+1,:)=temp;
        end
        if strcmp(inpfname,'solver_crack_ct_4_path1.inp')
            if ismember(i,[1,2,3,4])
                temp=pathElConn(i,:);
                temp=[temp(end) temp(1:end-1)];
                temp=[temp(end) temp(1:end-1)];
                temp=[temp(end) temp(1:end-1)];
                pathElConn(i,:)=temp;
            end
            if ismember(i+1,[43,44,45,46])
                temp=pathElConn(i+1,:);
                temp=[temp(end) temp(1:end-1)];
                pathElConn(i+1,:)=temp;
            end
        end
    end

    % Visualize Path
    pathCood=[];
    x=dataGauss.Pt;
    [sp_1, ~]=Q4ShapeFn(x(2,1),x(2,2));
    [sp_2, ~]=Q4ShapeFn(x(3,1),x(3,2));
    for i =1:size(pathElConn,1)
        el_connect=pathElConn(i,:);
        el_nodes=nodes(el_connect,:);
        x1=[sp_1(:)'*el_nodes(:,1) sp_1(:)'*el_nodes(:,2)];
        x2=[sp_2(:)'*el_nodes(:,1) sp_2(:)'*el_nodes(:,2)];
        gp=[x1;x2];
        plot(gp(:,1),gp(:,2));
        pathCood=[pathCood;x1;x2];
    end
%     plot(pathCood(:,1),pathCood(:,2), 'linewidth',1)
%     hold on

%     hold off
end