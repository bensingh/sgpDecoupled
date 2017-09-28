%Created by Bensingh Dhas  and Md Masiur Rahaman
%Descroption: Computes path integral for the path give using element number
%and the loacl edge number.
%Dependencies: None
%Created on: 23rd Nov, 2016
%Remark: Pls check the integration
%carefully.

function visoPlasJ= nearTipJIntegralAbaqusDynamic(pathElem,domainElem,connect,pathconnect,points,...
    disp_hat,disp_hatDot,disp_hat_dDot,comp_disp_hat,comp_disp_hatDot,comp_disp_hat_dDot,...
    gamma_p_hat,gamma_p_hat_dot,gamma_p_hat_dDot,omega_hat,omega_hat_dot,omega_hat_dDot,...
    displacement0,DisPdot0,gamma_p0,gamma_p_dot0,dataGauss,mat)

integrationPath= edgeNaturalUnitNormal(pathElem);

% Integrate along the path
wt=[1,1]; % since 1D quadrature.
jPath= [0;0];
path_len=0;
for i=1:size(integrationPath,1)
    node_index=pathconnect(i,1:4);
    el_nodes=points(node_index,:);
    
    unit_normal=integrationPath(i,6:7);
    el_gauss_pt=[integrationPath(i,2:3); integrationPath(i,4:5)];
    
    el_disp0=displacement0(node_index,:);
    el_gamma0=gamma_p0(node_index);
    el_gamma_dot0=gamma_p_dot0(node_index);
    
    el_disp_hat= disp_hat(node_index,:);
    el_comp_disp_hat=comp_disp_hat(node_index,:);
    
    el_gamma_hat=gamma_p_hat(node_index);
    el_gamma_hat_dot=gamma_p_hat_dot(node_index);
    el_omega_hat=omega_hat(node_index);
    el_omega_hat_dot=omega_hat_dot(node_index);
    
    
    [elJPath,el_path_len] = elementViscoplasticJIntegral_path(el_nodes,...
        el_disp_hat,el_comp_disp_hat,el_gamma_hat_dot,...
        el_omega_hat_dot,el_gamma_hat,el_omega_hat,el_disp0,el_gamma0,el_gamma_dot0,...
        el_gauss_pt,wt,unit_normal,mat);
    path_len=path_len+el_path_len;
    jPath=jPath+elJPath;
    
end

% Integrate over domain
jDomain= [0;0];

for i=1:numel(domainElem)
    it_el=domainElem(i);
    node_index=connect(it_el,1:4);
    el_nodes=points(node_index,:);
    
    el_disp_hat= disp_hat(node_index,:);
    el_disp_hatDot= disp_hatDot(node_index,:);
    el_disp_hat_dDot= disp_hat_dDot(node_index,:);
    
    el_comp_disp_hat=comp_disp_hat(node_index,:);
    el_comp_disp_hatDot=comp_disp_hatDot(node_index,:);
    el_comp_disp_hat_dDot=comp_disp_hat_dDot(node_index,:);
    
    el_DisPdot0 = DisPdot0(node_index,:);
    el_gamma0=gamma_p0(node_index,:);
    el_gamma_dot0=gamma_p_dot0(node_index);
    
    el_gamma_hat=gamma_p_hat(node_index);
    el_gamma_hat_dot=gamma_p_hat_dot(node_index);
    el_gamma_hat_dDot=gamma_p_hat_dDot(node_index);
    
    el_omega_hat=omega_hat(node_index);
    el_omega_hat_dot=omega_hat_dot(node_index);
    el_omega_hat_dDot=omega_hat_dDot(node_index);
    
    
     elJDomain =elementViscoplasticJIntegral_domainDynamic(el_nodes,el_disp_hat,el_disp_hatDot,el_disp_hat_dDot,...
        el_comp_disp_hat,el_comp_disp_hatDot,el_comp_disp_hat_dDot,el_gamma_hat,el_gamma_hat_dot,el_gamma_hat_dDot,...
        el_omega_hat,el_omega_hat_dot,el_omega_hat_dDot,el_DisPdot0,el_gamma0,el_gamma_dot0,dataGauss,mat);
    
    jDomain = jDomain + elJDomain;
    
end

visoPlasJ = jPath - jDomain;

end


