clear all;

addpath("Preprocess/");
addpath("Assembly/");
addpath("Solver/");
addpath("Modes_Plots/")
addpath("Meshes/");
addpath("Postprocess/");
addpath("AuxFunctions/");
addpath("SVD_Debug_Results/");


% Input data file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NAME_INPUT_DATA = 'DATA_ASSIGNMENT5' ;
%-------------------------------------------------

Ux=0.8;
mu=1.225e-3; %Dynamic viscosty [kg/(m s)]
rho=1.225; %Air density [kg/m3]
nu = mu/rho; %Kinematic Viscosity coefficient [m2/s] 
NavierStokes = 0.5;


% Definition of meshes 
% Name mesh
NameMeshP= 'LidDriven75'; %Cylinder10,20,40,75 LidDriven75 %NACA0012_AoA_5
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

% Load SVD Modes
load('SVD Matrices.mat')
Phi = SNAP_cluster.u.U;
Psi = SNAP_cluster.p.U;

n_mod = size(Phi,2);
Namecase=['_mu_',num2str(mu),'_U_',num2str(Ux)];

NameMeshV=strcat(NameMeshV,'NSS',Namecase);
NameMeshP=strcat(NameMeshP,'NSS',Namecase);

direc=['GIDPOST/',NameMeshP,'/'];
mkdir(direc);
NameMeshV=[direc,NameMeshV];
NameMeshP=[direc,NameMeshP];

u = Phi(:,1);
p = Psi(:,1);

ux = u;

u = ux(1:2:end-1);
v = ux(2:2:end);

GidPostProcess2DV(COOR_v,CN_v,TypeElement_v,u,v,NAME_INPUT_DATA,NameMeshV); 
GidPostProcess2DP(COOR_p,CN_p,TypeElement_p,p,NAME_INPUT_DATA,NameMeshP); 
