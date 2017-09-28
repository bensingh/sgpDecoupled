%Created by Bensingh Dhas
%Descroption: Reads abaqus .inp file and generates the nodes and connectivity table.
%Dependencies: None
%Created on: 18th, Oct, 2016.

function mesh =readMeshAbaqus(file_name)
    fid = fopen(file_name);
    tline = fgetl(fid);
    i=1;
    %Read the text stream into a cell line by li
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

    % Convert string data to number data
    % Read Node data
    node_data=readTextData(st_node+1,end_node-1,data,'table');
    mesh.nodes=node_data(:,2:3);
    mesh.noNd=size(node_data,1);

    % Read Connectivity data
    element_data=readTextData(st_element+1,end_element-1,data,'table');
    element_data=element_data(:,2:5);
    mesh.elements=element_data;
    mesh.noEl=size(element_data,1);

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