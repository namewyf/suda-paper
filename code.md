%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Example code for Localized Orthogonal Decomposition method in 2d
%  Supplementry material of the book
%% An Introduction to the Localized Orthogonal Decomposition method
%  by 
%% A. Malqvist (Chalmers University of Technology and University of Gothenburg)
%  and 
%% D. Peterseim (University of Augsburg)
%
% The code is explained in detail in Chapter 7 of the book. It produces an
% LOD approximation of a prototypical linear elliptic diffusion problem as
% used throughout the book. The script can be run in Matlab or Octave. 
% All required subroutines are contained in this file. The user may adapt 
% the parameters. 
%
% Copyright (C) 2020 Daniel Peterseim - All Rights Reserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Main script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clear work space
clear all;

%% Model problem
% physical dimension
d = 2; 
% seed the random number generator (if possible)
if exist('rng') rng(27); end 
% bounds of diffusion coefficient A
alpha = 0.01; 
beta = 1; 
% number of cells per coordinate direction
Neps = 2^6;
% diffusion coefficient on Cartesian mesh
Aeps = alpha+(beta-alpha).*(randi(2,Neps,Neps)-1);
% turn Aeps into function depending on (a list of) points in Omega
A = @(x) Aeps((ceil(x(:,1)*Neps)-1)*Neps+(ceil(x(:,2)*Neps)));
%A = @(x) 1;
% right-hand side
f = @(x) ones(size(x,1),1);

%% Initial mesh
% list of coordinates of vertices
p = [0 0; 1 0; 1 1; 0 1];
% triangles by references to vertex numbers
t = [1 2 4; 2 3 4];
% initial triangulation of Omega as structure array
T0 = struct('t',t,'p',p,'Dnodes',[1;2;3;4]);

%% LOD parameters
Hlevel =2;
hlevel = 8;
ell = 1;

%% coarse mesh
% refine initial mesh Hlevel times
[TH] = refineMesh(T0,Hlevel);
% define numbers of vertices NH and elements NTH
NH = size(TH.p,1);
NTH = size(TH.t,1);
% define number of local vertices
NHdg = numel(TH.t);
% enumerate local vertices acoording to their appearence in TH.t
dgHidx = reshape(1:NHdg,d+1,[]).';
% compute area of elements
areaTH = computeArea(TH);
% nngH is the number of neighbors for interior and 0 for boundary vertices
nngH = accumarray(TH.t(:),1,[NH,1]);
nngH(TH.Dnodes) = 0;

%% fine mesh
% refine initial mesh hlevel-Hlevel times
[Th,P1,R1,P0,P1dg] = refineMesh(TH,hlevel-Hlevel);
% define numbers of vertices Nh and elements NTh
Nh = size(Th.p,1);
NTh = size(Th.t,1);
% define number of local vertices
Nhdg = numel(Th.t);
% enumerate local vertices acoording to their appearence in Th.t
dghidx = reshape((1:Nhdg),d+1,[])';
% nngh is the number of neighbors for interior and 0 for boundary vertices
nngh = accumarray(Th.t(:),1,[Nh,1]);
nngh(Th.Dnodes) = 0;

%% quasi-interpolation
% embedding of continuous into discontinuous finite element functions
cg2dgh = sparse(1:Nhdg,reshape(Th.t.',1,[]),1,Nhdg,Nh);
% mass matrix for finite discontinuous finite elements
Mhdg = kron(spdiags(computeArea(Th),0,NTh,NTh),[2 1 1;1 2 1;1 1 2]./12);
% inverse of mass matrix for coarse discontinuous finite elements
B = kron(spdiags(1./computeArea(TH),0,NTH,NTH),[9 -3 -3;-3 9 -3;-3 -3 9]);
% compute L2 projection onto coarse discontinuous finite elements
PiHdg = B*((P1dg.'*Mhdg)*cg2dgh);
% embedding of continuous into discontinuous finite element functions
cg2dgH = sparse(1:NHdg,reshape(TH.t.',1,[]),1,NHdg,NH);
EH = spdiags(1./sum(cg2dgH,1).',0,NH,NH)*cg2dgH.';
% IH is the product of EH and PiHdg.
IH = EH*PiHdg;
%% for test
MHdg = kron(spdiags(computeArea(TH),0,NTH,NTH),[2 1 1;1 2 1;1 1 2]./12);

%% element patches
% incidence matrix of coarse elements and vertices
Ivt = sparse(TH.t,repmat((1:NTH).',1,d+1),1,NH,NTH);
% incidence matrix of coarse elements
Itt = spones(Ivt.'*Ivt);
% initialize patch matrix as identity
patch = speye(NTH);
% succesively apply Itt ell times
for k = 1:ell
    patch = Itt*patch;
end
% set non-zero entries to 1 and diagonal entries to 2
patch = spones(patch)+speye(NTH);

%% global stiffness and mass matrix
% evaluate diffusion coefficient at midpoints of elements
Ah = A((Th.p(Th.t(:,1),:)+Th.p(Th.t(:,2),:)+Th.p(Th.t(:,3),:))./3);
% assemble stiffness matrix for discontinuous functions
Shdg = assembleStiffnessMatrixDG(Th,Ah);

%% element correctors
% initialize cell array for element correctors
CT = cell(1,NTH);
% loop over coarse elements
for k=1:NTH
	% tpH(j) = true iff coarse element j is in ell-patch of element k 
	tpH = logical(patch(:,k));
	% indices of corase patch vertices
	dofpH = find(accumarray(reshape(TH.t(tpH,:),[],1),1,[NH 1]));
	% remove (global) boundary vertices	
	dofpH = dofpH(nngH(dofpH)>0); 
	% tph(j) = true iff fine element j is in ell-patch of element k 
	tph = logical(P0*patch(:,k));
	% nngph is the number of neighbors for fine patch vertices and 0 else
	nngph = accumarray(reshape(Th.t(tph,:),[],1),1,[Nh 1]);
	% non-zero entries of nngph identify fine patch vertices
	dofph = find(nngph);
	% remove boundary vertices by comparing number of neighbors 
	% in the patch with number of global neighbors
	dofph = dofph(nngph(dofph)==nngh(dofph));
	% restrict fine stiffness matrix to dof
	dofphdg = reshape(dghidx(tph,:).',[],1);
	cg2dghk = cg2dgh(dofphdg,dofph);
	Sph = (cg2dghk.'*Shdg(dofphdg,dofphdg))*cg2dghk;
    % tTh(j) = true iff fine element j is in coarse element k
    tTh = logical(P0*patch(:,k)>1);
    % assemble right-hand side
    dofThdg = reshape(dghidx(tTh,:).',[],1);
    cg2dghT = cg2dgh(dofThdg,dofph);
    rhsp = (cg2dghT.'*Shdg(dofThdg,dofThdg))*P1dg(dofThdg,dgHidx(k,:));
	% restrict quasi interpolation to patch
	IHp = IH(dofpH,dofph);
	% solve corrector problem
	X = Sph\[IHp.' rhsp];
	mu = (IHp*X(:,1:length(dofpH)))\(IHp*X(:,length(dofpH)+(1:d+1)));
     CT{k} = sparse(Nh,d+1);
     CT{k}(dofph,:) = X(:,length(dofpH)+(1:d+1))-X(:,1:length(dofpH))*mu;
     %% for test
%      mH=length(dofpH);
%      mh=length(dofph);
%      MHdgT=MHdg(dofpH,dofpH);
%      X=[Sph (MHdgT*IHp).';
%         MHdgT*IHp  zeros(mH,mH)];
%       u=X\[rhsp;zeros(mH,3)];
%  
% 
%       X=[Sph -(MHdgT*IHp).';
%         MHdgT*IHp  zeros(mH,mH)];
%       Tran=cg2dgH(:,dofpH);
%    Tran=Tran';
%    dof4Ek=3*k-2:3*k;
%    rhs2=Tran*MHdg(:,dof4Ek);
%       rhs=[zeros(mh,3);rhs2];
%       u=X\rhs;
%       CT{k}(dofph,:)=u(1:mh,:);
%      CTT{k} = sparse(Nh,d+1);
%      CTT{k}(dofph,:) = u(1:mh,:);
     
end
% global correction operator
C_ell = cell2mat(CT)*cg2dgH;
% LOD basis
G = P1-C_ell;
%G=cell2mat(CT)*cg2dgH;

%% Build and solve coarse LOD system
% coarse degrees of freedom (interior vertices)
dofH = true(NH,1);
dofH(TH.Dnodes) = 0;
% LOD stiffness matrix
% assemble stiffness and mass matrix for continuous finite elements
G0 = G(:,dofH);
SHLOD0 = ((G0.'*(cg2dgh.'*Shdg))*cg2dgh)*G0;
% right-hand side
rhs = G0.'*(cg2dgh.'*(Mhdg*(cg2dgh*(P1*f(TH.p)))));
% solve coarse linear system
uH = zeros(NH,1);
uH(dofH) = SHLOD0\rhs;
% reconstruct fine-scale information to compute LOD approximation
uHms = G*uH;

%% Reference solution
dofh = true(Nh,1);
dofh(Th.Dnodes) = 0;
cg2dgh0 = cg2dgh(:,dofh);
Sh = (cg2dgh.'*Shdg)*cg2dgh;
Mh = (cg2dgh.'*Mhdg)*cg2dgh;
uh = zeros(Nh,1);
uh(dofh) = Sh(dofh,dofh)\(Mh(dofh,:)*f(Th.p));

%% Vizualization
% We shall finally plot the reference solution along with the LOD
% approximation. 

%% test
% figure(1)
% huHh = trisurf(Th.t,Th.p(:,1),Th.p(:,2),uHms,'EdgeColor','none','FaceColor','interp');
% camproj perspective;campos([-3,-5,1]);daspect([1,1,0.5]);box on;set(gca,'BoxStyle','full');
% light;lighting phong;camlight('left');
% axis([0 1 0 1 0 0.5]);caxis([0 .3]);

figure
set(gcf, 'Position',  [10, 10, 1200, 300])
% plot coefficient
subplot(1,4,1)
[X,Y] = meshgrid(linspace(0,1,Neps),linspace(0,1,Neps));
hA = surf(X,Y,Aeps,'EdgeColor','none');
hold on; axis equal; view(0,90);colormap(gca,gray);caxis([alpha,beta]);
h = title({'Diffusion coefficient $A$',['($\alpha = ',num2str(alpha),'$ (black) , $\beta = ',num2str(beta),'$ (white))']});
set(h,'Interpreter','latex')
% plot reference solution
subplot(1,4,2)
huh = trisurf(Th.t,Th.p(:,1),Th.p(:,2),uh,'EdgeColor','None');
camproj perspective;campos([-3,-5,1]);daspect([1,1,0.5]);box on;set(gca,'BoxStyle','full');
light;lighting phong;camlight('left');
axis([0 1 0 1 0 0.5]);caxis([0 .3]);
h = title({['Reference solution $u_h$'],['($h=2^{-',num2str(hlevel-0.5),'}$)']});
set(h,'Interpreter','latex')
% plot LOD approximation
subplot(1,4,3)
huHh = trisurf(Th.t,Th.p(:,1),Th.p(:,2),uHms,'EdgeColor','none','FaceColor','interp');
camproj perspective;campos([-3,-5,1]);daspect([1,1,0.5]);box on;set(gca,'BoxStyle','full');
light;lighting phong;camlight('left');
axis([0 1 0 1 0 0.5]);caxis([0 .3]);
h = title({'Full LOD approximation $u^{ms,h}_{H,\ell}$',['($h=2^{-',num2str(hlevel-0.5),'}$, $H=2^{-',num2str(Hlevel-0.5),'}$, $\ell=',num2str(ell),'$)']});
set(h,'Interpreter','latex')
% plot FE part of LOD approximation
subplot(1,4,4)
 huHh = trisurf(TH.t,TH.p(:,1),TH.p(:,2),uH,'EdgeColor','k','LineWidth',0.5,'FaceColor','interp');
camproj perspective;campos([-3,-5,1]);daspect([1,1,0.5]);box on;set(gca,'BoxStyle','full');
light;lighting phong;camlight('left');
axis([0 1 0 1 0 0.5]);caxis([0 .3]);
h = title({['FE-part of LOD approximation $u^{h}_{H,\ell}$'],['($h=2^{-',num2str(hlevel-0.5),'}$, $H=2^{-',num2str(Hlevel-0.5),'}$, $\ell=',num2str(ell),'$)']});
set(h,'Interpreter','latex')
% % plot difference of LOD approximation and reference solution
% subplot(1,4,4)
% huHh = trisurf(Th.t,Th.p(:,1),Th.p(:,2),uh-uHms,'EdgeColor','none','FaceColor','interp');
% camproj perspective;campos([-3,-5,1]);daspect([1,1,0.5]);box on;set(gca,'BoxStyle','full');
% light;lighting phong;camlight('left');


%% Error
% Error of the LOD approximation in the energy norm
fprintf('LOD-Error (in energy norm): %7.5e\n',sqrt((uh-uHms).'*(Sh*(uh-uHms))));
% Error of the LOD approximation in the L^2 norm
fprintf('LOD-Error (in L2 norm): %7.5e\n',sqrt((uh-uHms).'*(Mh*(uh-uHms))));
% Error of FE-part og the LOD approximation in the L^2 norm
fprintf('Error of FE-part of LOD (in L2 norm): %7.5e\n',sqrt((uh-P1*uH).'*(Mh*(uh-P1*uH))));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Subroutines 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Function refineMesh 
% refines a coarse triangulation TH nref times by applying the 
% function REFINE nref times. Suitable restriction and prolongation
% operator between corresponding mesh functions are computed
function [Th,P1,R1,P0,P1dg] = refineMesh(Th,nref)
    % number of vertices
    NH = size(Th.p,1);
    % compute fine reference mesh
    P1 = speye(NH);
    P0 = speye(size(Th.t,1));
    P1dg = speye(numel(Th.t));
    for k=1:nref
        [Th,p,p0,p1dg] = refine(Th);
        P1dg = p1dg*P1dg;
        P1 = p*P1;
        P0 = p0*P0;
    end
    % compute coarse mesh nodal interpolant
    Nh = size(Th.p,1);
    [i,j] = find(P1>1-2*eps);
    R1 = sparse(j,i,ones(length(i),1),NH,Nh);
end
% End of function refineMesh

%% Function refine
% refines the current triangulation by dividing
% each marked element into 2^dim congruent elements (red refinement)
function [T,P,P0,P1dg] = refine(T)
    % Construct data structure
    [np] = size(T.p,1); 
    nt = size(T.t,1);
    [e,nmbe] = computeEdges(T);
    ne = size(e,1);
    d2p = sparse(e,e(:,[2,1]),(1:ne)'*[1 1],np,np);
    % New nodes from the mid points of each edge
    newnode = 0.5.*(T.p(e(:,1),:)+T.p(e(:,2),:)); 
    P = sparse((1:ne)'*[1 1],e,.5,ne,np);
    T.p = [T.p; newnode]; 
    P = [speye(np);P];
    emarker = (np+1:np+ne)';
    p = [T.t,emarker(d2p(T.t(:,1)+np*(T.t(:,2)-1))),...
             emarker(d2p(T.t(:,2)+np*(T.t(:,3)-1))),...
             emarker(d2p(T.t(:,3)+np*(T.t(:,1)-1)))];
    T.t = [p(:,[1,4,6]);p(:,[4,2,5]);p(:,[6,5,3]);p(:,[4,5,6])];
    tnew2t = repmat(1:nt,1,4).';
    E = speye(nt);     
    P1dg = [kron(E,[1 0 0;.5 .5 0;.5 0 .5]);kron(E,[.5 .5 0;0 1 0;0 .5 .5]);
            kron(E,[.5 0 .5;0 .5 .5;0 0 1]);kron(E,[.5 .5 0;0 .5 .5;.5 0 .5])];
    Dnodes = false(np,1);
    Dnodes(T.Dnodes) = true;
    T.Dnodes = find([Dnodes;Dnodes(e(:,1)) & Dnodes(e(:,2)) & nmbe==1]);
    P0 = sparse((1:length(tnew2t))',tnew2t,1,length(tnew2t),nt);
end
% End of function REFINE

%% Function computeEdges 
% computes the edges of the triangulation T
function [e,nmbe] = computeEdges(T)
    np = size(T.p,1);
    e = [T.t(:,[1,2]); T.t(:,[1,3]); T.t(:,[2,3])];
    e = sort(e,2);
    d2p = sparse(e(:,1),e(:,2),1,np,np);
    [e1,e2,nmbe] = find(d2p);
    e = [e1,e2];
end
% End of function computeEdges

%% Function assembleStiffnessMatrixDG calculates stiffness matrix with 
% respect to simplicial mesh Th with optional Th-piecewise constant scalar 
% diffusion coefficient coeff
% discontinuous Galerkin
function [S] = assembleStiffnessMatrixDG(Th,coeff)
    % compute gradients of FE functions
    ve1 = Th.p(Th.t(:,3),:)-Th.p(Th.t(:,2),:);
    ve2 = Th.p(Th.t(:,1),:)-Th.p(Th.t(:,3),:);
    ve3 = Th.p(Th.t(:,2),:)-Th.p(Th.t(:,1),:);
    area2 = abs(ve3(:,2).*ve2(:,1)-ve3(:,1).*ve2(:,2));
    G(:,:,1)=[-ve1(:,2) ve1(:,1)]./(area2*[1 1]);
    G(:,:,2)=[-ve2(:,2) ve2(:,1)]./(area2*[1 1]);
    G(:,:,3)=[-ve3(:,2) ve3(:,1)]./(area2*[1 1]);
    c = coeff.*computeArea(Th);
    % local stiffness matrices
    S = zeros(size(Th.t,1),9);
    for j = 1:3
        for k = 1:3
            S(:,k+3*(j-1)) = c.*dot(G(:,:,j),G(:,:,k),2);
        end
    end
    % create global matrix
    I = (1:3).'*ones(1,3);
    ndof = numel(Th.t);
    idx = reshape(1:ndof,3,[]).';
    S = sparse(idx(:,reshape(I.',1,9)),idx(:,reshape(I,1,9)),S,ndof,ndof);
end
% End of function assembleStiffnessMatrix

%% Function computeArea
% calculates the area of the elements of the triangulation T
function v = computeArea(T)
    d12 = T.p(T.t(:,2),:)-T.p(T.t(:,1),:);
    d13 = T.p(T.t(:,3),:)-T.p(T.t(:,1),:);
    v = (d12(:,1).*d13(:,2)-d12(:,2).*d13(:,1))./2;
end
% End of function computeArea