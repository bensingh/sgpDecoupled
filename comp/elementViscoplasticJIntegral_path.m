
%Authors: Md Masiur Rahaman and Bensingh Dhas 
%Description: Solves the gurtin mode of visco-plasticity.
%Load control with micro inertia.
%Microforce balance as the outer loop.
%Created on: 29 th Aug, 2017

function [el_incr_J_int_path,path_len] = elementViscoplasticJIntegral_path(el_nodes,...
    el_disp_hat,el_comp_disp_hat,el_gamma_hat_dot,el_omega_hat_dot,...
    el_gamma_hat,el_omega_hat,el_disp, el_gamma,el_gamma_dot,el_gauss_pt,gauss_wt,unit_normal,mat)

m1=mat.m1;
l1=mat.l1;
D = mat.D;
h0=mat.H; s0 = mat.S_0;
nu=mat.nu;
d_0 = mat.d_0 ;

[sp,~]=Q4ShapeFn(el_gauss_pt(1,1)/abs(el_gauss_pt(1,1)),el_gauss_pt(1,2)/abs(el_gauss_pt(1,2)));
gp1=[sp*el_nodes(:,1) sp*el_nodes(:,2)];
[sp,~]=Q4ShapeFn(el_gauss_pt(2,1)/abs(el_gauss_pt(2,1)),el_gauss_pt(2,2)/abs(el_gauss_pt(2,2)));
gp2=[sp*el_nodes(:,1) sp*el_nodes(:,2)];
path_len=norm(gp1-gp2);


el_incr_J_int_path =[0;0];

for it_gp=1:size(el_gauss_pt,1)
    
    x=el_gauss_pt(it_gp,:);
    B= Q4StrainDisplacement(x(1),x(2),el_nodes);
    incr_strain_vector=B*el_disp_hat(:);
    incr_stress_vector=D*incr_strain_vector;
    [~,dsp]=Q4ShapeFn(x(1),x(2));
    jacobian=dsp*el_nodes;
    dsp = inv(jacobian)*dsp;
    
    gamma_hat_dot=sp*el_gamma_hat_dot(:);
    omega_hat_dot=sp*el_omega_hat_dot(:);
    gamma_hat=sp*el_gamma_hat(:);
    omega_hat=sp*el_omega_hat(:);
    
    gamma_dot=sp*el_gamma_dot(:);
    gamma=sp*el_gamma(:);
    grad_gamma=dsp*el_gamma(:);
    
    gamma_dot(gamma_dot<=1e-10) = 1e-10;
    
    c1=s0*((h0/s0)*(gamma_dot/d_0)^m1);
    c2 =s0*((1+(h0/s0)*gamma)*(m1/d_0^m1)*gamma_dot^(m1-1));
    pi_0 =s0*(1+(h0/s0)*gamma)*(gamma_dot/d_0)^m1;

    
    du_dx(1,1)=dsp(1,:)*el_disp_hat(:,1);
    du_dx(1,2)=dsp(2,:)*el_disp_hat(:,1);
    du_dx(2,1)=dsp(1,:)*el_disp_hat(:,2);
    du_dx(2,2)=dsp(2,:)*el_disp_hat(:,2);

    du_comp_dx(1,1)=dsp(1,:)*el_comp_disp_hat(:,1);
    du_comp_dx(1,2)=dsp(2,:)*el_comp_disp_hat(:,1);
    du_comp_dx(2,1)=dsp(1,:)*el_comp_disp_hat(:,2);
    du_comp_dx(2,2)=dsp(2,:)*el_comp_disp_hat(:,2);
    
    
    strain=B*el_disp(:);
    stress=D*strain;
    sig_3D = [stress(1)  stress(3) 0;  stress(3) stress(2) 0 ; 0  0  nu*(stress(1)+stress(2))];
    inc_sig_3D = [incr_stress_vector(1)  incr_stress_vector(3) 0;  incr_stress_vector(3) incr_stress_vector(2) 0 ;...
        0  0  nu*(incr_stress_vector(1)+incr_stress_vector(2))];
    sig_3D_dev = sig_3D - (trace(sig_3D)/3)*eye(3,3);
    inc_sig_3D_dev = inc_sig_3D - (trace(inc_sig_3D)/3)*eye(3,3);
    tau = sqrt(sum(sum(sig_3D_dev.^2)));
    stress_dev_vec=[sig_3D_dev(1,1), sig_3D_dev(2,2), sig_3D_dev(1,2)];
    inc_stress_dev_vec=[inc_sig_3D_dev(1,1), inc_sig_3D_dev(2,2), inc_sig_3D_dev(1,2)];
    grad_gamma_hat=dsp*el_gamma_hat(:);
    grad_omega_hat=dsp*el_omega_hat(:);
    gp_omega_hat=sp*el_omega_hat(:);
    
    incr_stress_tensor=[incr_stress_vector(1) incr_stress_vector(3);incr_stress_vector(3) incr_stress_vector(2)];
    incr_traction=incr_stress_tensor*unit_normal(:);
    %%-----------------------------first term-----------------------%%
    
    incr_lag=0.5*(incr_strain_vector'*D*incr_strain_vector) + ...
        + s0*l1^2* grad_gamma_hat'*grad_omega_hat + c1*gamma_hat* omega_hat ...
        + 0.5*c2*(gamma_hat_dot*omega_hat-omega_hat_dot*gamma_hat)-...
        inc_stress_dev_vec*stress_dev_vec'*(omega_hat/tau);
    
    %% terms due to last known state
    incr_lag=incr_lag+(stress'*incr_strain_vector)+(s0*l1^2*grad_gamma(:)'*grad_omega_hat(:))-(tau-pi_0)*omega_hat;
    
    first_term = incr_lag*unit_normal(:);
    
    %%-----------------------------second term-----------------------%%
    
    stress_bar_vec=(D*stress_dev_vec(:))/tau;
    stress_bar_tensor=[stress_bar_vec(1) stress_bar_vec(3);stress_bar_vec(3) stress_bar_vec(2)];
    traction_bar=stress_bar_tensor*unit_normal(:);

    micro_traction = (grad_gamma_hat*(grad_omega_hat'* unit_normal(:)) + grad_omega_hat*(grad_gamma_hat'* unit_normal(:)));
    
    second_term = 0.5*(du_dx'+du_comp_dx')*incr_traction + s0*l1^2*micro_traction - gp_omega_hat*du_dx'*traction_bar;
    
     %% terms due to last known state 
    stress_mat = [stress(1)  stress(3) ;  stress(3) stress(2)]; 
    traction_0 = stress_mat*unit_normal(:);
    incr_strain_mat = [incr_strain_vector(1) incr_strain_vector(3); incr_strain_vector(3) incr_strain_vector(2)];
     
    second_term = second_term + (incr_strain_mat*traction_0) + (s0*l1^2*grad_gamma*grad_gamma_hat'*unit_normal(:));
    el_incr_J_int_path =el_incr_J_int_path +gauss_wt(it_gp)*path_len*(first_term - second_term);
    
end

end