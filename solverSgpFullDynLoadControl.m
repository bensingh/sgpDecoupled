%Authors: Md Masiur Rahaman and Bensingh Dhas
%Description: Solves the gurtin mode of simplified visco-plasticity.
%Load control with micro and macro inertia.
%Created on: 31Aug, 2017

close all;
clear all; clc;
addpath('./bin');
addpath('./comp');
% Read Complete data from abaqus inp
inpfname='./inp/solver_crack_ct_4_path1_symm.inp';
[nodes,connect,unit_force_vector,dead_dof_u,dead_dof_gp,pathJ,domainJ,track_node]=readNodesConnectAbaqus(inpfname);
mesh.nodes=nodes;
mesh.elements=connect;
mesh.noEl=size(connect,1);
mesh.noNd=size(nodes,1);

dataGauss=gaussData(2);
meshGp = deformedGausspt(mesh,dataGauss);
xp=squeeze(meshGp(:,:,1));
yp=squeeze(meshGp(:,:,2));
% Material model
mat=materialModelDynamic(1);

% Problem information
noNd=size(nodes,1);
noEl=size(connect,1);
string=strcat('No Nodes: ', num2str(noNd), 'No Elements: ',num2str(noEl));
disp(string);

%Input data used to geneterate plots
crackLineNodes=find(mesh.nodes(:,2)==0);
crackTip=112;
xAxisIndex=find(abs(mesh.nodes(:,2))<0.00001);
[~,xIndex]=sort(mesh.nodes(xAxisIndex,1));
xAxisIndexSorted=xAxisIndex(xIndex(:,1));

trackNode=[29 7];

% Forcing data
t=0; tmax= 2; n_t = 5e2; dt= tmax/n_t;
load_max=5e6; %=>N/m^2
rate_loading = load_max/tmax;

% Nodal states
displacement=zeros(noNd,2);
velocity=zeros(noNd,2);
accleration=zeros(noNd,2);

gammaP=zeros(noNd, 1);
gammaPdot=zeros(noNd, 1);
gammaPdDot=zeros(noNd, 1);
convGammaP=zeros(noNd, 1);

% State at gauss point
noGp=2;
gpStress=zeros(noEl, noGp^2, 3);
trialGpStress=zeros(noEl, noGp^2, 3);
gpEpsilonP0=zeros(noEl, noGp^2, 3);
gpNodes=zeros(noEl, noGp^2, 2);

% Initialize data containers
dataDisplacement={};
dataDispDot={};
dataDisPdDot={};

dataGammaPhat={};
datagammaPhatDot={};
datagammaPhat_dDot={};

dataGammaP={};
dataGammaPdot={};
dataGammaPdDot={};

dataGPStress={};

% Elastic Stiffness and the deviatoric elastic stiffness
kElas=spalloc(2*noNd,2*noNd,8*noNd);

for itEl=1:noEl
  nodeIndex=connect(itEl,1:4);
  elNodes=nodes(nodeIndex,:);
  elDof=[nodeIndex nodeIndex+noNd];
  keE=Q4PlaneDynaElasticOperator(elNodes,dataGauss,mat,dt);
  kElas(elDof,elDof)=kElas(elDof,elDof)+keE;
end


activeDof=setdiff(1:noNd*2,dead_dof_u(:));
itTime=1;

while t<tmax % Time loop
    t=t+dt;
    msg=strcat('Time step: ',num2str(t), '\n');
    fprintf(msg);
    
    noLoadDof=length(find(unit_force_vector>0));
    fTotal=(t*rate_loading/noLoadDof)*unit_force_vector(1:2*noNd);

    gpEpsilonP = gpEpsilonP0;
    
    tsDisp = displacement;
    tsVel = velocity;
    tsAccler = accleration;
    
    
    %% Elastic prediction

    fInternal=zeros(2*noNd,1);
    
    for itEl=1:noEl
        
        nodeIndex=connect(itEl,1:4);
        elNodes=nodes(nodeIndex,:);
        elDof=[nodeIndex nodeIndex+noNd];
        elDisp=tsDisp(nodeIndex,:);
        elAcc=tsAccler(nodeIndex,:);
        elEpsilonP=squeeze(gpEpsilonP(itEl,:,:));
        fe=sgpDynaInternalForce(elNodes,elDisp,elAcc,elEpsilonP,dataGauss,mat);
        fInternal(elDof,1)=fInternal(elDof,1)+fe;
        
    end
    
    elasResidue=fTotal-fInternal;    
    errElas = 100;        


    while errElas>1e-4 
        
        sol=zeros(2*noNd,1);
        sol(activeDof,1)=kElas(activeDof,activeDof)\elasResidue(activeDof,1);
        tsDisp=tsDisp+ reshape(sol,noNd,2);
        tsAccler = (4/dt^2)*(tsDisp-(displacement + dt*velocity))-accleration;
        tsVel = velocity + (dt/2)*(tsAccler + accleration);

        % Compute Stress
        for itEl=1:mesh.noEl
            nodeIndex=mesh.elements(itEl,1:4);
            elDof=[nodeIndex nodeIndex+mesh.noNd];
            elNodes=mesh.nodes(nodeIndex,:);
            elDisp=tsDisp(nodeIndex,:);
            elEpsilonP=squeeze(gpEpsilonP(itEl,:,:));
            gpStress(itEl,:,:)= sgpElementStress(elNodes,elDisp,elEpsilonP,dataGauss,mat);
        end
        
            
        itgammaP=gammaP;
        itgammaPdot=gammaPdot;
        itgammaPdDot=gammaPdDot;
        
        intMicroForce=residueDynMicoForce(mesh,itgammaP,itgammaPdot,itgammaPdDot,gpStress,dataGauss,mat);
        errGamma=100;
        itNr=1;
        
        while errGamma>1e-4 
            % Construct Tangent stiffness for micro force

            ke=planeDynPlasOperator(mesh,itgammaP,itgammaPdot,dt,dataGauss,mat);
            sol=(ke\-intMicroForce);
            itgammaP= itgammaP+sol;
            itgammaPdDot = (4/dt^2)*(itgammaP-(gammaP + dt*gammaPdot))-gammaPdDot;
            itgammaPdot = gammaPdot + (dt/2)*(itgammaPdDot + gammaPdDot);
            itgammaP(itgammaP<=0) = 0;
            itgammaPdot(itgammaPdot<=0) = 0;
            intMicroForce=residueDynMicoForce(mesh,itgammaP,itgammaPdot,itgammaPdDot,gpStress,dataGauss,mat);
            errGamma=norm(intMicroForce);
            msg=strcat('     Residue Norm: ',num2str(errGamma));
            disp(msg);
            itNr=itNr+1;

        end
  
        gammaPhat = itgammaP - gammaP;
        gammaPhatDot = itgammaPdot - gammaPdot;
        gammaPhat_dDot = itgammaPdDot - gammaPdDot;
        
        gammaP=itgammaP;
        gammaPdot=itgammaPdot;
        gammaPdDot=itgammaPdDot;

 
        % Update plastic strain
        for itEl=1:noEl
            nodeIndex=connect(itEl,1:4);
            elDof=[nodeIndex nodeIndex+noNd];
            elNodes=nodes(nodeIndex,:);
            elDisp0=displacement(nodeIndex,:);
            elGammaDot0= gammaPdot(nodeIndex);
            elGammaDot=itgammaPdot(nodeIndex);
            elEpsilonP0=squeeze(gpEpsilonP0(itEl,:,:));
            gpEpsilonP(itEl,:,:)=sgpPlasticStrain(elNodes,elDisp0,elGammaDot0,elGammaDot,elEpsilonP0,dt,dataGauss,mat);
        end
        
        % Elastic residue
        fInternal=zeros(2*noNd,1);
        
        for itEl=1:noEl
            
            nodeIndex=connect(itEl,1:4);
            elNodes=nodes(nodeIndex,:);
            elDof=[nodeIndex nodeIndex+noNd];
            elDisp=tsDisp(nodeIndex,:);
            elAcc=tsAccler(nodeIndex,:);
            elEpsilonP=squeeze(gpEpsilonP0(itEl,:,:));  %% approximation because of inexact tangent
            fe=sgpDynaInternalForce(elNodes,elDisp,elAcc,elEpsilonP,dataGauss,mat);
            fInternal(elDof,1)=fInternal(elDof,1)+fe;
            
        end
        
        elasResidue=fTotal-fInternal;        
        errElas=norm(elasResidue(activeDof));
        msg=strcat('  Elast Residue Norm: ',num2str(errElas));
        disp(msg);

    end
    
    
    displacement=tsDisp;
    velocity=tsVel;
    accleration=tsAccler;
    
    gpEpsilonP0 = gpEpsilonP;


if mod(itTime,10)==0
    
        zp1=squeeze(gpStress(:,:,1));
        zp2=squeeze(gpStress(:,:,2));
        zp3=squeeze(gpStress(:,:,3));
        figure(1)
        scatter(xp(:),yp(:),25,zp1(:),'filled');
        colorbar
        figure(2)
        scatter(xp(:),yp(:),25,zp2(:),'filled');
        colorbar
        figure(3)
        scatter(xp(:),yp(:),25,zp3(:),'filled');
        colorbar
    figure(4)
    scatter(nodes(:,1),nodes(:,2),25,gammaP,'filled');
    colorbar
%     caxis([0 1e-6]);
    pause (0.001)
    
end

  dataDisplacement{itTime}=displacement;
  dataDispDot{itTime}=velocity;
  dataDisPdDot{itTime}=accleration;
  
  
  dataGammaPhat{itTime}=gammaPhat;
  datagammaPhatDot{itTime}=gammaPhatDot;
  datagammaPhat_dDot{itTime}=gammaPhat_dDot;
  
  dataGammaP{itTime}=gammaP;
  dataGammaPdot{itTime}=gammaPdot;
  dataGammaPdDot{itTime}=gammaPdDot;
  
  
  dataGPStress{itTime}=gpStress;

  p_delta(itTime,1)=t;
  p_delta(itTime,2)=displacement(trackNode(1),2)-displacement(trackNode(2),2);
  p_delta(itTime,3)=gammaP(crackTip);
  
  itTime=itTime+1;
  
  msg=strcat('Time: ', num2str(t), ' Time Step: ', num2str(itTime), '\n');

end



%% ---------------Complementary solution---------------
noTimeit=itTime-1
active_op_dof=1:noNd;

% Solution of complementary equation

[dataOmegaPhat,dataOmegaPhat_dot,dataOmegaPhat_dDot]= complementaryDynamic(nodes,connect,...
  active_op_dof,p_delta,dataGammaPhat,datagammaPhatDot,datagammaPhat_dDot,dataGammaP,dataGammaPdot,dataGauss,dt,mat);


%% --------------Calculation of J-Integral---------------
jIntegral=zeros(noTimeit-1,4);
jIntegral(:,1)=linspace(0,t-dt,noTimeit-1);
for i=1:3
  inpfname=strcat('./inp/solver_crack_ct_4_path',num2str(i),'.inp');
  [nodes,connect,unit_force_vector,dead_dof_u,dead_dof_gp,pathJ,domainJ,track_node]=readNodesConnectAbaqus(inpfname);
  JpathElConn =plotJintegral_path(inpfname, dataGauss);
  
  for it_time=1:noTimeit-1
    disp(it_time)
    
    disp_hat = dataDisplacement{it_time+1}-dataDisplacement{it_time};
    disp_hatDot = dataDispDot{it_time+1}-dataDispDot{it_time};
    disp_hat_dDot = dataDisPdDot{it_time+1}-dataDisPdDot{it_time};
    
    comp_disp_hat = disp_hat;
    comp_disp_hatDot = -disp_hatDot;
    comp_disp_hat_dDot = disp_hat_dDot;
    
    displacement0=dataDisplacement{it_time};
    DisPdot0=dataDispDot{it_time};
    
    gamma_p0 = dataGammaP{it_time};
    gamma_p_dot0 = dataGammaPdot{it_time};
    
    gamma_p_hat = dataGammaP{it_time+1}-dataGammaP{it_time};
    gamma_p_hat_dot=dataGammaPdot{it_time+1}-dataGammaPdot{it_time};
    gamma_p_hat_dDot=dataGammaPdDot{it_time+1}-dataGammaPdDot{it_time};
    
    % Pls think and correct it
    omega_hat=dataOmegaPhat{end-it_time};
    omega_hat_dot=dataOmegaPhat_dot{end-it_time};
    omega_hat_dDot=dataOmegaPhat_dDot{end-it_time};
    
    inc_j_path= nearTipJIntegralAbaqusDynamic(pathJ,domainJ,connect,JpathElConn,nodes,...
      disp_hat,disp_hatDot,disp_hat_dDot,comp_disp_hat,comp_disp_hatDot,comp_disp_hat_dDot,...
      gamma_p_hat,gamma_p_hat_dot,gamma_p_hat_dDot,omega_hat,omega_hat_dot,omega_hat_dDot,...
      displacement0,DisPdot0,gamma_p0,gamma_p_dot0,dataGauss,mat);
    
    jIntegral(it_time+1,i+1)=inc_j_path(1);
    
  end
  
end

gpdat=[mesh.nodes(xAxisIndexSorted,1) gammaP(xAxisIndexSorted)];
gpXfname=strcat('in_gpXaxis_m_0.10','.csv');
csvwrite(gpXfname, gpdat);

gpfname=strcat('in_gammap_m_0.10','.csv');
csvwrite(gpfname, gammaP);

gpHfname=strcat('in_gp_Hist_m_0.10','.csv');
csvwrite(gpHfname, p_delta(:,2:3));

jfname=strcat('in_jintegralData_m_0.10','.csv');
csvwrite(jfname, jIntegral);

figure(6)
hold on
plot(p_delta(:,2)*1e3,jIntegral(:,2),'-*');
plot(p_delta(:,2)*1e3,jIntegral(:,3),'-*');
plot(p_delta(:,2)*1e3,jIntegral(:,4),'-*');
hold off
%% ------------------------END of the programme------------------------%%
