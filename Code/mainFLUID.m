%%% Finite Element Solver for the Navier-Stokes Flow  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------
% Technical University of Catalonia
% Charbel Mouawad and Ferran 
% Conducted by Prof. Dr. Joaquin Hernández Ortega
% -------------------------------------------------------------------------

% Start the profiler
profile on

clear
close all

addpath("Preprocess/");
addpath("Assembly/");
addpath("Solver/")
addpath("Meshes/")
addpath("Postprocess/");
addpath("AuxFunctions/");
addpath(genpath("ML"));
%% INPUT  %% 
% Input data file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NAME_INPUT_DATA = 'DATA_ASSIGNMENT5' ;
%-------------------------------------------------


%% PREPROCESS 
% Inputs
% ---------
mu=0.1; %Dynamic viscosty [kg/(m s)]
rho=1.225; %Air density [kg/m3]
nu = mu/rho; %Kinematic Viscosity coefficient [m2/s] 
Ux=1; %Selects Boundary conditions for the kovasznay problem
NavierStokes=1; %Select 1 for NS, 0 for Stokes, 0.5 for stepped NS, 1.5 for Stepped Stokes
res=1e-9; %Select residual for NS
maxiter=1e+3; %Maximum iterations for NS, if not converged exit with error
rel_factor=0.8; %Relaxation factor, <1 for under relaxation >1 for over relaxation !!!VERY IMPORTANT IF STEPPED!!! IF NO CONVERGENCE IS REACHED, TRY TO REDUCE MESH SIZE
debug=0; %Debug parameter for enabling extra functions. 0 disabled, 1 enabled.
steps=10; %Number of stepping for solving nonlinear problem
tolU=1e-6; %SVD tolerance for pressure and velocity
tolP=1e-6; %SVD tolerance for pressure and velocity
nModes = 7;


% ----------

% Definition of meshes 
% Name mesh
NameMeshP= 'Cylinder40'; %Cylinder10,20,40,75 LidDriven75 %NACA0012_AoA_5
NameMeshV=[NameMeshP '_v'];
NameMeshNodes=[NameMeshP, '_MODES_'];
% Velocity mesh
NameFileMeshDATA = ['./Meshes/' NameMeshV];
[COOR_v,CN_v,TypeElement_v,TypeElementB_v, ViscMglo,rnod_v,dR_v,...  
    tracglo_v,CNb_v,~,~] = ReadInputDataFile_v(NameFileMeshDATA,Ux,nu); 
% Pressure mesh
NameFileMeshDATA = ['./Meshes/' NameMeshP];
[COOR_p,CN_p,TypeElement_p,TypeElementB_p, ~,rnod_p,dR_p,...  
    tracglo_p,CNb_p,~,~] = ReadInputDataFile_p(NameFileMeshDATA); 
% ----------


%%Names for Post processing
Namecase=['_mu_',num2str(mu),'_U_',num2str(Ux),'_Res_',num2str(res),'_Steps_',num2str(steps),'_tolU_',num2str(tolU),'_tolP_',num2str(tolP)];
switch(NavierStokes)
    case(1)
    NameMeshV=strcat(NameMeshV,'NS',Namecase);
    NameMeshP=strcat(NameMeshP,'NS',Namecase);
    case(0)
    NameMeshV=strcat(NameMeshV,'S',Namecase);
    NameMeshP=strcat(NameMeshP,'S',Namecase);
    case(0.5)
    NameMeshV=strcat(NameMeshV,'NSS',Namecase);
    NameMeshP=strcat(NameMeshP,'NSS',Namecase);
    case(1.5)
    NameMeshV=strcat(NameMeshV,'SS',Namecase);
    NameMeshP=strcat(NameMeshP,'SS',Namecase);
end
direc=['GIDPOST/',NameMeshP,'/'];
mkdir(direc);
NameMeshV=[direc,NameMeshV];
NameMeshP=[direc,NameMeshP];

% Number of dimensions is obtained from the mesh
nDOFTOT=size(COOR_p,1)+size(COOR_v,1)*2;
disp(['Total number of DOFs: ',num2str(nDOFTOT)])

totalTime=tic;

%% SOLVER
% Computation of K and G matrices
disp('Computing KGL')
[K,G,F,Bst,OmegaGlo] = ComputeKGF(COOR_v,CN_v,COOR_p,CN_p,CNb_v,tracglo_v,TypeElement_v, TypeElement_p,TypeElementB_v,nu,ViscMglo,debug);


% Solution of the system of equations
switch(NavierStokes)
    case(1)
    disp('Solving the system of equations for Navier Stokes all at once')
    [u,v,p] = SolverNavierStokes(COOR_v,CN_v,rnod_v,dR_v,COOR_p,rnod_p,dR_p,K,G,res,TypeElement_v,maxiter,rel_factor,Bst,OmegaGlo,debug,nModes);
    case(0.5)
     disp('Solving the system of equations for Navier Stokes in steps')
    % [u,v,p,SNAP_cluster,SNAP] = SolverNavierStokesSteps(COOR_v,CN_v,rnod_v,dR_v,COOR_p,rnod_p,...
    %     dR_p,K,G,res,TypeElement_v,maxiter,rel_factor,Bst,OmegaGlo,debug,steps,NameMeshP,tolU,tolP);
    [ux,v,p,SNAP_cluster,SNAP] = SVDtotal_2(COOR_v,CN_v,rnod_v,dR_v,COOR_p,rnod_p,dR_p,K,G,res,...
    TypeElement_v,maxiter,rel_factor,Bst,OmegaGlo,debug,steps,NameMeshP,tolU,tolP);
    case(0)
    disp('Solving the system of equations for Stokes')
    [u,v,p,ul,pl] = SolverStokes(COOR_v,rnod_v,dR_v,COOR_p,rnod_p,dR_p,K,G,F); 
    case (1.5)
    disp('Solving the system of equations of Stokes in Steps')
    [SNAP_cluster,SNAP] = SolverStokesSteps(COOR_v,CN_v,rnod_v,dR_v,COOR_p,rnod_p,dR_p,K,G,res,...
    TypeElement_v,maxiter,rel_factor,Bst,OmegaGlo,debug,steps,NameMeshP,tolU,tolP,F);   
end
totalTime=toc(totalTime);
disp(['Total time to solve the problem ', num2str(totalTime),' s']);
%% POST-PROCESS
disp('Starting the postprocessing')

% addpath("SVD_Debug_Results/");
% load("SVD Matrices.mat");


if (NavierStokes == 1 || NavierStokes == 0)
    GidPostProcess2DV(COOR_v,CN_v,TypeElement_v,u,v,NAME_INPUT_DATA,NameMeshV); 
    GidPostProcess2DP(COOR_p,CN_p,TypeElement_p,p,NAME_INPUT_DATA,NameMeshP); 

elseif (NavierStokes == 1.5 || NavierStokes == 0.5)

    Phi = SNAP_cluster.u.U;
    Psi = SNAP_cluster.p.U;

    sol = Phi(:,1);
    p = Psi(:,1);

    u = sol(1:2:end-1);
    v = sol(2:2:end);
    u = u';
    v = v';

    GidPostProcess2DV(COOR_v,CN_v,TypeElement_v,u,v,NAME_INPUT_DATA,NameMeshV); 
    GidPostProcess2DP(COOR_p,CN_p,TypeElement_p,p,NAME_INPUT_DATA,NameMeshP); 

end

%Postproc(COOR_v,COOR_p,u,v,p) %FOR LID DRIVEN ONLY

% Stop the profiler
profile off

% View the profiling result
profile viewer