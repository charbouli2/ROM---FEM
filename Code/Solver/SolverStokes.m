function [u,v,p,ul,pl] = SolverStokes(COOR_v,rnod_v,dR_v,COOR_p,rnod_p,dR_p,K,G,F)
% Solves the linear system of equations that result from Stokes flow
%Inputs:
%   - COOR_v: Matrix with the coordiantes of velocity nodes
%   - COOR_p: Matrix with the coordiantes of pressure nodes
%   - rnod_v: Vector containing the restricted nodes of the velocity
%   - rnod_p: Vector containing the restricted nodes of the pressure
%   - dR_v: Vector containing the prescribed velocities
%   - dR_p: Vector containing the prescribed pressures
%   - K: Global viscosity matrix 
%   - G: Global gradient operator
%Outputs:
%   - u: vector with the horizontal component of the velocity
%   - v: vector with the vertical component of the velocity
%   - p: vector with the pressure
    %  ("svd_stokes.mat");
    
    load("svd_stokes_Cylinder.mat")
    Phi = SNAP_cluster.u.U;
    Psi = SNAP_cluster.p.U;
    rank_p = size(Psi,2);

    ndim = size(COOR_v,2); 
    % Application of Boundary Conditions
    nDOF_v = ndim * size(COOR_v,1);
    DOFl_v = (1:nDOF_v)';
    DOFr_v = rnod_v;
    DOFl_v(DOFr_v) = [];
    
    % Pressure on node 1 set to zero
    nDOF_p = size(COOR_p,1);
    DOFl_p = (1:nDOF_p)';
    DOFr_p = rnod_p;
    DOFl_p(DOFr_p) = [];
    
    
    % Decomposition of matrices
    Kll = K(DOFl_v,DOFl_v);
    Gl = G(DOFl_p,DOFl_v);
    GlT = Gl';
    Gr = G(DOFl_p,DOFr_v);
    Klr = K(DOFl_v,DOFr_v);
    L = zeros(rank_p);
    % Fl=F(DOFl_v);
    
    Phi = [Phi, GlT*Psi];

    A=[Phi'*Kll*Phi Phi'*GlT*Psi; Psi'*Gl*Phi L];
    B=[-Phi'*Klr*dR_v; -Psi'*Gr*dR_v];
    
    % Solution of the system
    d = A\B;
    % Decomposition of the different terms

    q = d(1:size(Phi,2));
    e = d(size(Phi,2)+1:end);
    % 
    % A1 = Phi'*GlT*Psi;
    % B1 = -Phi'*Klr*dR_v - Phi'*Kll*Phi*q;
    % 
    % e1 = A1\B1;

    ul = Phi*q;
    pl = Psi*e;

    sol(DOFl_v) = ul;
    sol(DOFr_v) = dR_v;
    
    p=zeros(size(COOR_p,1),1);
    p(DOFl_p)=pl;
    p(DOFr_p)=dR_p;
    u = sol(1:2:end-1);
    v = sol(2:2:end);
    
end