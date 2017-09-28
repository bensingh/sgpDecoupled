
function mat = materialModelDynamic(Material)

switch Material
    
    case 1   %% C(T) specimen
        
        mat.rho = 0; %%7850;
        mat.nu   = 0.3;
        mat.E    = 200e9;
        mat.G    = mat.E/(2*(1+mat.nu));     % Shear Modulus
        mat.K    = 2*mat.G*(1+mat.nu)/(3*(1-2*mat.nu));  % Bulk Modulus
        mat.lam  = (2* mat.G*mat.nu)/(1-2*mat.nu);
        mat.S_0  = 270e6;  % Initial Yield Stress in mata
        mat.d_0  = 1e-4;% Reference Strain Rate
        mat.m1    = 0.1;  % Strain rate sensitivity parameter
        mat.l1   = 1e-5;    % Energetic length scale
        mat.H    = 500e6;   % Hardening Modulus in mata
        mat.D = (mat.E/((1+mat.nu)*(1-2*mat.nu)))*[1-mat.nu mat.nu 0 ; mat.nu 1-mat.nu 0; 0 0 (1-2*mat.nu)/2];
        mat.lm = sqrt(1e-5);
     case 2   %% C(T) specimen
        
        mat.rho = 0; %%7850;
        mat.nu   = 0.3;
        mat.E    = 200e9;
        mat.G    = mat.E/(2*(1+mat.nu));     % Shear Modulus
        mat.K    = 2*mat.G*(1+mat.nu)/(3*(1-2*mat.nu));  % Bulk Modulus
        mat.lam  = (2* mat.G*mat.nu)/(1-2*mat.nu);
        mat.S_0  = 270e6;  % Initial Yield Stress in mata
        mat.d_0  = 1e-4;% Reference Strain Rate
        mat.m1    = 0.12;  % Strain rate sensitivity parameter
        mat.l1   = 1e-5;    % Energetic length scale
        mat.H    = 500e6;   % Hardening Modulus in mata
        mat.D = (mat.E/((1+mat.nu)*(1-2*mat.nu)))*[1-mat.nu mat.nu 0 ; mat.nu 1-mat.nu 0; 0 0 (1-2*mat.nu)/2];
        mat.lm = sqrt(1e-5);
      case 3   %% C(T) specimen
        
        mat.rho = 0; %%7850;
        mat.nu   = 0.3;
        mat.E    = 200e9;
        mat.G    = mat.E/(2*(1+mat.nu));     % Shear Modulus
        mat.K    = 2*mat.G*(1+mat.nu)/(3*(1-2*mat.nu));  % Bulk Modulus
        mat.lam  = (2* mat.G*mat.nu)/(1-2*mat.nu);
        mat.S_0  = 270e6;  % Initial Yield Stress in mata
        mat.d_0  = 1e-4;% Reference Strain Rate
        mat.m1    = 0.15;  % Strain rate sensitivity parameter
        mat.l1   = 1e-5;    % Energetic length scale
        mat.H    = 500e6;   % Hardening Modulus in mata
        mat.D = (mat.E/((1+mat.nu)*(1-2*mat.nu)))*[1-mat.nu mat.nu 0 ; mat.nu 1-mat.nu 0; 0 0 (1-2*mat.nu)/2];
        mat.lm = sqrt(1e-5);
       
    otherwise
        
        warning('Material_model : ! Unexpected Material Type !');
        
end

end
