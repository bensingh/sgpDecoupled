%Created by Bensingh Dhas
%Descroption: Reads abaqus .inp file and generates the nodes and connectivity table.
%Dependencies: None
%Created on: 18th, Oct, 2016.

function [node_data,element_data,unit_force_vector,dead_dof_u,dead_dof_gp,pathJ,domainJ,track_node] =readNodesConnectAbaqus(file_name)
    dof_pn=3;
    fid = fopen(file_name);
    tline = fgetl(fid);
    i=1;
    %Read the text stream into a cell line by line.
    while ischar(tline)
        data{i,1}=tline;
        tline = fgetl(fid);
        i=i+1;
    end
    %Find the line number where each data ends.
    st_node =find(ismember(data,'Node'));
    end_node =find(ismember(data,'EndNode'));
    st_element =find(ismember(data,'Element'));
    end_element=find(ismember(data,'EndElement'));
    st_load =find(ismember(data,'ForceBc'));
    end_load =find(ismember(data,'EndForceBc'));
    st_disp =find(ismember(data,'DispBc'));
    end_disp =find(ismember(data,'EndDispBc'));
    st_pathJ =find(ismember(data,'PathJ'));
    end_pathJ=find(ismember(data,'EndPathJ'));
    st_domainJ =find(ismember(data,'DomainJ'));
    end_domainJ=find(ismember(data,'EndDomainJ'));
    st_track_node =find(ismember(data,'TrackNode'));
    end_track_node=find(ismember(data,'EndTrackNode'));

    % Convert string data to number data
    % Read Node data
    node_data=readTextData(st_node+1,end_node-1,data,'table');
    node_data=node_data(:,2:3);
    no_nodes=size(node_data,1);

    % Read Connectivity data
    element_data=readTextData(st_element+1,end_element-1,data,'table');
    element_data=element_data(:,2:5);

    % Create index of dead dof

    %     bottomFixedNodes=find(node_data(:,2)<-0.299);
    %     dead_dof_u=[bottomFixedNodes(:); no_nodes+bottomFixedNodes(:)];
    %     dead_dof_gp=2*no_nodes+bottomFixedNodes(:);

    % Create force vector
    %     topLoadNodes=find(node_data(:,2)>0.29990);
    %     unit_force_vector=zeros(dof_pn*no_nodes,1);
    %     force_dof_index=no_nodes+topLoadNodes;
    %     unit_force_vector(force_dof_index)=1;

    % Track node
    %     track_node=topLoadNodes;
    % Read load data
    load_data=readTextData(st_load+1,end_load-1,data,'table');
    unit_force_vector=zeros(dof_pn*no_nodes,1);
    force_dof_index=no_nodes*(load_data(:,2)-1)+load_data(:,1);
    unit_force_vector(force_dof_index)=load_data(:,3);

    % Read Displacement BC Data
    disp_bc_data=readTextData(st_disp+1,end_disp-1,data,'table');
    dead_x=disp_bc_data(find(disp_bc_data(:,2)==1),1);
    dead_y=disp_bc_data(find(disp_bc_data(:,3)==1),1)+no_nodes;
    dead_gp=disp_bc_data(find(disp_bc_data(:,4)==1),1)+2*no_nodes;
    dead_dof_u=[dead_x(:);dead_y(:)];
    dead_dof_gp=dead_gp(:);

    % Read Path for J integral calculation
    pathJ=readTextData(st_pathJ+1,end_pathJ-1,data,'table');
    
    % Read Domain for J integral calculation
    domainJ=readTextData(st_domainJ+1,end_domainJ-1,data,'stream');

    % Read nodes to be tracked
    track_node=readTextData(st_track_node+1,end_track_node-1,data,'stream');


end

function section_data=readTextData(start_line,end_line,data,option)
    j=1;
    for i=start_line:end_line
        temp=data{i,1};
        section_data{j,1}=temp(find(~isspace(temp))); %wicked code think rewriting
        section_data{j,1}=str2num(section_data{j,1});
        j=j+1;
    end 
    if strcmp(option,'table')
        section_data=cell2mat(section_data);
    end
    if strcmp(option,'stream')
        section_data=[section_data{:}];
    end
end