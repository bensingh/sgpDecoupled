

function el_incr_J_int_domain =elementViscoplasticJIntegral_domainDynamic(el_nodes,el_disp_hat,el_disp_hatDot,el_disp_hat_dDot,...
    el_comp_disp_hat,el_comp_disp_hatDot,el_comp_disp_hat_dDot,el_gamma_hat,el_gamma_hat_dot,el_gamma_hat_dDot,...
    el_omega_hat,el_omega_hat_dot,el_omega_hat_dDot,el_DisPdot0,el_gamma0,el_gamma_dot0,dataGauss,mat)

s0=mat.S_0; 
h0=mat.H; d_0=mat.d_0;
m1=mat.m1;
lm = mat.lm ; rho = mat.rho;

    el_incr_J_int_domain =[0;0];

    for it_gp=1:numel(dataGauss.Wt)

        x=dataGauss.Pt(it_gp,:);
        [sp,dsp]=Q4ShapeFn(x(1),x(2));
        jacobian=dsp*el_nodes;
        dsp = inv(jacobian)*dsp;
        JxW=dataGauss.Pt(it_gp)*det(jacobian);
        
        u_hat_dot=[sp*el_disp_hatDot(:,1);sp*el_disp_hatDot(:,2)];
        u_hat_dDot=[sp*el_disp_hat_dDot(:,1); sp*el_disp_hat_dDot(:,2)];
        v_hat_dot= [sp*el_comp_disp_hatDot(:,1);sp*el_comp_disp_hatDot(:,2)];
        v_hat_dDot=[sp*el_comp_disp_hat_dDot(:,1);sp*el_comp_disp_hat_dDot(:,2)];
      
        gamma_hat=sp*el_gamma_hat(:);
        gamma_hat_dot=sp*el_gamma_hat_dot(:);
        gamma_hat_dDot=sp*el_gamma_hat_dDot(:);
         
        omega_hat=sp*el_omega_hat(:);
        omega_hat_dot=sp*el_omega_hat_dot(:);
        omega_hat_dDot=sp*el_omega_hat_dDot(:);
        
        grad_u_hat=dsp*el_disp_hat;
        grad_u_hat_dot=dsp*el_disp_hatDot;
        grad_v_hat=dsp*el_comp_disp_hat;
        grad_v_hat_dot=dsp*el_comp_disp_hatDot;

        grad_gamma_hat=dsp*el_gamma_hat(:);
        grad_gamma_hat_dot=dsp*el_gamma_hat_dot(:);
        
        grad_omega_hat=dsp*el_omega_hat(:);
        grad_omega_hat_dot=dsp*el_omega_hat_dot(:);

        DisPdot = [sp*el_DisPdot0(:,1);sp*el_DisPdot0(:,2)];
        gamma=sp*el_gamma0(:);
        gammaPdot=sp*el_gamma_dot0(:);
        gammaPdot(gammaPdot<=1e-10) = 1e-10;

        c2 =s0*((1+(h0/s0)*gamma)*(m1/d_0^m1)*gammaPdot^(m1-1));
        c1_dyn = rho*DisPdot;
        c2_dyn = rho*lm^2*gammaPdot;
        
        el_incr_J_int_domain =el_incr_J_int_domain +JxW*0.5*c2*...
            (grad_omega_hat_dot*gamma_hat + grad_omega_hat * gamma_hat_dot...
            -grad_gamma_hat_dot*omega_hat - grad_gamma_hat* omega_hat_dot)...
             -JxW*grad_v_hat_dot*c1_dyn - JxW*rho*...
            (grad_v_hat*u_hat_dDot+grad_v_hat_dot*u_hat_dot...
            +grad_u_hat*v_hat_dDot+grad_u_hat_dot*v_hat_dot)...
            -JxW*c2_dyn*grad_omega_hat_dot - JxW*rho*lm^2*...
            (gamma_hat_dDot*grad_omega_hat+gamma_hat_dot*grad_omega_hat_dot...
            +omega_hat_dDot*grad_gamma_hat_dot+omega_hat_dot*grad_gamma_hat_dot);

    end

end
