function CH = fluidHomogenization(lx, ly, zeta, phi, x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lx        = Unit cell length in x-direction.
% ly        = Unit cell length in y-direction.
% zeta      = material parameters. Two entries.
% phi       = Angle between horizontal and vertical cell wall. Degrees
% x         = Material indicator matrix. Size used to determine nelx/nely
% Example: homogenize(1,2,[2.1e5,1.2e5],[0.2,0.4],30,randi(2,10,10));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INITIALIZE
% Deduce discretization
[nely, nelx] = size(x);
% Stiffness matrix consists of two parts, one belonging to lambda and
% one belonging to mu. Same goes for load vector
dx = lx/nelx; dy = ly/nely;
nel = nelx*nely;
[ke, ke_brink,be,pe,le] = elementMatVec(dx/2, dy/2,phi);
% Node numbers and element degrees of freedom for full (not periodic) mesh
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nel,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nel,1);
edofVecp = 2*(nelx+1)*(nely+1) + reshape(1*nodenrs(1:end-1,1:end-1)+1,nel,1);
edofMatp = repmat(edofVecp,1,4) + repmat([0 nely+[1 0] -1],nel,1);
%% IMPOSE PERIODIC BOUNDARY CONDITIONS
% Use original edofMat to index into list with the periodic dofs
nn = (nelx+1)*(nely+1); % Total number of nodes
nnP = (nelx)*(nely);    % Total number of unique nodes
nnPArray = reshape(1:nnP, nely, nelx);
% Extend with a mirror of the top border
nnPArray(end+1,:) = nnPArray(1,:);
% Extend with a mirror of the left border
nnPArray(:,end+1) = nnPArray(:,1);
% Make a vector into which we can index using edofMat:
dofVector = zeros(2*nn, 1);
dofVector(1:2:end) = 2*nnPArray(:)-1;
dofVector(2:2:end) = 2*nnPArray(:);
edofMat = dofVector(edofMat);
dofVector(2*nn+1:3*nn) = 2*nnP+nnPArray(:);
edofMatp= dofVector(edofMatp);
ndof = 3*nnP; % Number of dofs
%% ASSEMBLE STIFFNESS MATRIX
% Indexing vectors
zeta = zeta(1)*(x==0) + zeta(2)*(x==1);
K = repmat(ke(:),1,nel)+ke_brink(:)*zeta(:).';
B = repmat(be(:),nel,1);
P = repmat(pe(:),nel,1);
iK = kron(edofMat,ones(8,1))';
iB = iK(:,1:2:end);
iP = kron(edofMatp,ones(4,1))';
jK = kron(edofMat,ones(1,8))';
jB = kron(edofMatp,ones(1,8))';
jP = jB(1:2:end,:);
sA = [K(:); B; B; -P];
iA = [iK(:); iB(:); jB(:); iP(:)];
jA = [jK(:); jB(:); iB(:); jP(:)];
A = sparse(iA,jA,sA);
%% LOAD VECTORS AND SOLUTION
% Assembly three load cases corresponding to the three strain cases
iF = [reshape(edofMat(:,1:2:end)',4*nel,1)' reshape(edofMat(:,2:2:end)',4*nel,1)']';
jF = [ones(4*nel,1); 2*ones(4*nel,1)];
sF = repmat(le(:),2*nel,1);
F = sparse(iF,jF,sF,3*nnP,2);
% Solve (remember to constrain one node)
chi = zeros(ndof,2);
solfor = 1:ndof-1; % all except last pressure dof
chi(solfor,:) = A(solfor,solfor)\F(solfor,:);
%% HOMOGENIZATION
% The displacement vectors corresponding to the unit strain cases
CH = zeros(2);
CH(1,1) = sum(chi(dofVector(1:2:2*nn),1));
CH(1,2) = sum(chi(dofVector(2:2:2*nn),1));
CH(2,1) = sum(chi(dofVector(1:2:2*nn),2));
CH(2,2) = sum(chi(dofVector(2:2:2*nn),2));
CH = CH/nel;
disp('--- Homogenized elasticity tensor ---'); disp(CH)

%% COMPUTE ELEMENT STIFFNESS MATRIX AND FORCE VECTOR (NUMERICALLY)
function [ke, ke_brink, be, pe, le] = elementMatVec (a,b,phi)
% Constitutive matrix contributions
CMu = diag([2 2 1]); CLambda = zeros(3); CLambda(1:2,1:2) = 1; 
% Two Gauss points in both directions
xx=[-1/sqrt(3), 1/sqrt(3)]; yy = xx;
ww=[1,1];
% Initialize
ke = zeros(8); ke_brink = zeros(8);
be = zeros(8,4); pe = zeros(4);
le = zeros(4,1);
h2 = 4*(a^2+b^2+2*a*b*abs(cos(phi)/sin(phi)));
stab = h2/12;
L = zeros(3,4); L(1,1) = 1; L(2,4) = 1; L(3,2:3) = 1;
for ii=1:length(xx)
  for jj=1:length(yy)
    % Integration point
    x = xx(ii); y = yy(jj);
    N = 1/4*[(1-y)*(1-x) 0 (1-y)*(1+x) 0 (1+y)*(1+x) 0 (1-x)*(1+y) 0; ...
            0 (1-y)*(1-x) 0 (1-y)*(1+x) 0 (1+y)*(1+x) 0 (1-x)*(1+y)] ;
    % Differentiated shape functions
    dNx = 1/4*[-(1-y)  (1-y) (1+y) -(1+y)];
    dNy = 1/4*[-(1-x) -(1+x) (1+x)  (1-x)];
    % Jacobian
    J = [dNx; dNy]*[-a a a+2*b/tan(phi*pi/180) 2*b/tan(phi*pi/180)-a; ...
        -b -b b b]';
    detJ = J(1,1)*J(2,2) - J(1,2)*J(2,1);
    invJ = 1/detJ*[J(2,2) -J(1,2); -J(2,1) J(1,1)];
    Bp = invJ*[dNx;dNy; ];
    % Weight factor at this point
    weight = ww(ii)*ww(jj)*detJ;
    % Strain-displacement matrix
    G = [invJ zeros(2); zeros(2) invJ];
    dN = zeros(4,8);
    dN(1,1:2:8) = dNx;
    dN(2,1:2:8) = dNy;
    dN(3,2:2:8) = dNx;
    dN(4,2:2:8) = dNy;
    B = L*G*dN;
    % Element matrices
    ke = ke + weight*(B'*CMu*B);
    ke_brink = ke_brink + weight*(N'*N);
    be = be + weight*(B'*[1 1 0]'*N(1:4:end));
    pe = pe + weight*(Bp'*Bp*stab);
    le = le + weight*(N(1,1:2:end)');
  end
end