function [SNAP_cluster_tot,SNAP_tot] = SolverStokesSteps(COOR_v,CN_v,rnod_v,dR_v,COOR_p,rnod_p,dR_p,K,G,res,...
    TypeElement_v,maxiter,rel_factor,Bst,OmegaGlo,debug,steps,NameMeshP,tolU,tolP,F)
%Inputs:
%   - COOR_v: Matrix with the coordiantes of velocity nodes
%   - COOR_p: Matrix with the coordiantes of pressure nodes
%   - rnod_v: Vector containing the restricted nodes of the velocity
%   - rnod_p: Vector containing the restricted nodes of the pressure
%   - dR_v: Vector containing the prescribed velocities
%   - dR_p: Vector containing the prescribed pressures
%   - K: Global viscosity matrix 
%   - G: Global gradient operator
%   - res: tolerance of the residual for the nonlinear system of equations
%   - TypeElement_v: Velocity type of element for computing C
%   - maxiter: maximum iterations before the nonlinear solver returns an
%   error
%   - rel_factor: relaxation factor for the nonlinear solver
%   - Bst: B-stacked matrix for velocity
%   - OmegaGlo: Matrix containing the products of all Jacobians and weights
%   at all Gauss points
%   - debug: Debug parameter to control the assembly of C
%Outputs:
%   - u: vector with the horizontal component of the velocity
%   - v: vector with the vertical component of the velocity
%   - p: vector with the pressure
    ndim = size(COOR_v,2); 
    nnodeE_v = size(CN_v,2) ;

    nDOF_v = ndim * size(COOR_v,1);
    DOFl_v = (1:nDOF_v)';
    DOFr_v = rnod_v;
    DOFl_v(DOFr_v) = [];
    
    % Pressure on node 1 set to zero
    nDOF_p = size(COOR_p,1);
    DOFl_p = (1:nDOF_p)';
    DOFr_p = rnod_p;
    DOFl_p(DOFr_p) = [];
    
    %For PostProcessing
    TypeIntegrand = 'K';
    [~,posgp,~,~] = ComputeElementShapeFun(TypeElement_v,nnodeE_v,TypeIntegrand) ;
    
    % Decomposition of matrices
    Kll = K(DOFl_v,DOFl_v);
    Gl = G(DOFl_p,DOFl_v);
    GlT = Gl';
    Gr = G(DOFl_p,DOFr_v);
    Klr = K(DOFl_v,DOFr_v);
    %Auxiliary matrix
    L = zeros(length(DOFl_p));
    Fl=F(DOFl_v);

    
    %Initialization of SVD terms
    icluster=1;
    SNAP_tot.u=zeros(nDOF_v,steps); %% Edited : Changed the dimension (para solo calcular SVD de no velocidades no conocidas: Numero de grados de libertad dl)
    SNAP_tot.p=zeros(nDOF_p,steps); %% Edited : Changed the dimension (para solo calcular SVD de no presiones no conocidas: Numero de grados de libertad pl)

    SNAP.u=zeros(size(DOFl_v,1),steps); %% Edited : Changed the dimension (para solo calcular SVD de no velocidades no conocidas: Numero de grados de libertad dl)
    SNAP.p=zeros(size(DOFl_p,1),steps); %% Edited : Changed the dimension (para solo calcular SVD de no presiones no conocidas: Numero de grados de libertad pl)

    fff = fieldnames(SNAP_tot) ;
    for iii = 1:length(fff)
            SNAP_cluster.(fff{iii}) = [] ;
            SNAP_cluster.(fff{iii}).U = [] ;
            SNAP_cluster.(fff{iii}).S = [] ;
            SNAP_cluster.(fff{iii}).V = [] ;
    end
    
    VAR_tot.u=zeros(nDOF_v,1); %% Edited : Changed the dimension (para solo calcular SVD de no velocidades no conocidas: Numero de grados de libertad dl)
    VAR_tot.p=zeros(nDOF_p,1);    %% Edited : Changed the dimension (para solo calcular SVD de no presiones no conocidas: Numero de grados de libertad pl)
    
    VAR.u=zeros(size(DOFl_v,1),1); %% Edited : Changed the dimension (para solo calcular SVD de no velocidades no conocidas: Numero de grados de libertad dl)
    VAR.p=zeros(size(DOFl_p,1),1);    %% Edited : Changed the dimension (para solo calcular SVD de no presiones no conocidas: Numero de grados de libertad pl)

    DATA.STEPS=0:1/steps:1;
    DATA.STORE.NSTEPS_CLUSTER={1:1:steps};
    DATA.STORE.NAME_MATFILE_STORE={[NameMeshP,'.mat']};
    DATA.STORE.VAR.u=1;
    DATA.STORE.VAR.p=1;
    DATA.STORE.COMPRESS_WITH_SVD=1;
    DATA.STORE.TOLERANCE_SVD_COMPRESSION.u=tolU;
    DATA.STORE.TOLERANCE_SVD_COMPRESSION.p=tolP;

    %Solve of the nonlinear system
    disp('Solving Navier-Stokes in increments...')
    % Allocate ones for first iteration  
    deltadRv=dR_v/steps; %Increments
    dR_v=dR_v/steps; %Initial dRv
    
    istep=1; % $$ Hay que mirar si esto esta bien o tiene que ser istep = 1 $$
    time1=tic;
    while istep<=steps

        dR_v=dR_v+deltadRv;
        
        A=[Kll GlT; Gl L];
        B=[Fl-Klr*dR_v; -Gr*dR_v];

        d = A\B;

        ul = d(1:length(DOFl_v));
        pl = d(size(Kll,1)+1:end);

        
        CONVERGED=1;
        
        VAR_tot.u(DOFr_v)=dR_v;
        VAR_tot.u(DOFl_v)=ul;
        VAR_tot.p(DOFl_p)=pl; % Modifié by me, tu as mis pl a la place de p

        VAR.u = ul;
        VAR.p = pl;

        [SNAP_tot,DATA,icluster,SNAP_cluster_tot] = StoreInfoSnapshots(istep,icluster,SNAP_tot,VAR_tot,DATA,CONVERGED,COOR_v,CN_v,posgp,TypeElement_v,DOFl_v,NameMeshP,SNAP_cluster) ;
        
        [SNAP,DATA,icluster,SNAP_cluster] = StoreInfoSnapshots(istep,icluster,SNAP,VAR,DATA,CONVERGED,COOR_v,CN_v,posgp,TypeElement_v,DOFl_v,NameMeshP,SNAP_cluster) ;


        istep = istep + 1;
    end
             
    time1=toc(time1);
    disp(['Time to solve: ',num2str(time1),'s With ', num2str(steps),' steps'])
   
    save("NACA_Stokes_TOT.mat","SNAP_cluster_tot", "SNAP_tot");
    save("NACA_Stokes.mat","SNAP_cluster", "SNAP");
   
end