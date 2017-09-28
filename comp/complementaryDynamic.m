%Authors: Bensingh Dhas and Md Masiur Rahaman
%Description: Solves the complemetary equation and stores the data history.
%Dependencies: None
%Created on: 2nd Nov, 2016
function  [dataOmegaPhat,dataOmegaPhat_dot,dataOmegaPhat_dDot]= complementaryDynamic(nodes,connect,...
    active_dof,p_delta,dataGammaPhat,datagammaPhatDot,datagammaPhat_dDot,dataGammaP,dataGammaPdot,dataGauss,dt,mat)    

    dataOmegaPhat={};
    dataOmegaPhat_dot={};
    dataOmegaPhat_dDot={};
    
    no_nodes=size(nodes,1);
    no_elements=size(connect,1);
    
    % Saving the values dofs at the first time instant
    dataOmegaPhat{1}=dataGammaPhat{end};
    dataOmegaPhat_dot{1}=-datagammaPhatDot{end};
    dataOmegaPhat_dDot{1}= datagammaPhat_dDot{end}; %%zeros(numel(dataOmegaPhat{1}),1); 
    
    % Initialize the time Zero state of Complementary equation
    
    for it_t=1:size(p_delta,1)-1
          delta_t=p_delta(end-(it_t-1),1)-p_delta(end-it_t,1);
         msg=strcat('Time Step No. :  ',num2str(it_t) ,'; ','Delta t:',num2str(delta_t),'\n');
        fprintf(msg);
        
        gamma_p0=dataGammaP{end-(it_t-1)};
        gamma_p_dot0=dataGammaPdot{end-(it_t-1)};

        
        omegaP0=dataGammaPhat{end-(it_t-1)};
        omegaPdot0=-datagammaPhatDot{end-(it_t-1)};
        omegaPdDot0= datagammaPhat_dDot{end-(it_t-1)}; %%zeros(numel(omegaP0),1); 
          
        kg_comp=zeros(no_nodes,no_nodes);
        fg_comp=zeros(no_nodes,1);
        compSol=zeros(no_nodes,1);
        
        for it_el=1:no_elements
            
            % Get element data
            node_index=connect(it_el,1:4);
            el_nodes=nodes(node_index,:);
            el_gamma0=gamma_p0(node_index);  
            el_gamma_dot0=gamma_p_dot0(node_index);  
            
            el_omegaP0=omegaP0(node_index);  
            el_omegaPdot0=omegaPdot0(node_index);  
            el_omegaPdDot0=omegaPdDot0(node_index);
            [ke_comp,fe_comp]=elementComplementaryDynamic(el_nodes,el_gamma0,el_gamma_dot0,el_omegaP0,el_omegaPdot0,el_omegaPdDot0,dt,dataGauss,mat);
            
            % Assembly of global stiffness matrix
            kg_comp(node_index,node_index)=kg_comp(node_index,node_index)+ke_comp;
            fg_comp(node_index,1)=fg_comp(node_index,1)+fe_comp;

        end
      
        % Recover displacements and incremental Omega_p dot
        compSol(active_dof,1)=kg_comp(active_dof,active_dof)\ (fg_comp(active_dof));
        dataOmegaPhat{it_t+1}=omegaP0+compSol;
        dataOmegaPhat_dDot{it_t+1} = (4/dt^2)*(dataOmegaPhat{it_t+1}-(omegaP0 + dt*omegaPdot0))- omegaPdDot0;  %% check this
        dataOmegaPhat_dot{it_t+1} = omegaPdot0 + (2/dt)*(dataOmegaPhat{it_t+1}-(omegaP0 + dt*omegaPdot0));

    end
