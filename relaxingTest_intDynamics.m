% relaxingTest_intDynamics
% This script compares the hydrodynamic bead model against a relaxing rod
% probelm solved using a finite difference scheme.
% Finite Difference scheme is central divided differences in space, with
% implicit (backward) Euler in time.

% edited in Atom and uploaded using GitHub desktop.
% extra edit done in Atom editor.

clear all; close all
% viewInt = 0;

%% Parameter setup - finite difference method.

xMin = 0;
xMax = 1;
dx   = 2e-3;
x    = [xMin:dx:xMax];
Nx   = length(x);

q = 1/2;
gam = (2*log(2*q/0.005) + 1)^(1/4);

%% Parameter setup - bead model.

a = 0.005;                                       % radius of sphere approximating each bead; non-dimensional.
mu = 1e-3;                                      % dynamic viscosity of surrounding fluid.

L = 50e-6;                                      % length scale.
calB = 0;                                       % brownian dimensionless group (zero in experiment).
beta = 1;
calS = (1/(beta))*5e4;                        % spring-bending dimensionless group.

Nb71 = 71;                                        % number of beads.

epsilon = a;

tMin = 0;
tMax = 5e-5;

dt   = 1e-9;
t    = [tMin:dt:tMax];
Nt   = length(t);

Nt = length(t);                                  % number of time steps.
b071 = 1/(Nb71-1);

%% Set initial position.

% for finite difference method:
ySave      = zeros(Nx,2);                        % save current time step and one previous.
ySave(:,1) = intConfig(x);
yc(1)      = mean(ySave(:,1));

% for bead model:
xbSave71    = zeros(3,Nb71,2);                       % stores bead pos. at 2 time steps.
xc71        = zeros(3,Nt);                         % stores the xyz positions of the centre of mass at nt time-steps.

xbSave71(1,:,1) = linspace(0,L,Nb71)/L;              % equally space nodes in x-plane.
xbSave71(2,:,1) = intConfig(xbSave71(1,:,1));
xc71(:,1)       = mean(xbSave71(:,:,1),2);               % initial centre of mass.

%% Setup FD scheme matrix.

M   = zeros(Nx);
sig = gam^4*(dt/(dx^4));

R1  = horzcat((1+2*sig), -4*sig, +2*sig, zeros(1,Nx-3));
R2  = horzcat(-2*sig, (1+5*sig), -4*sig, +1*sig, zeros(1,Nx-4));
RN1 = horzcat(zeros(1,Nx-4), 1*sig, -4*sig, (1+5*sig), -2*sig);
RN  = horzcat(zeros(1,Nx-3), 2*sig, -4*sig, (1+2*sig));

M = diag((1+6*sig).*ones(1,Nx)) + diag(-4*sig.*ones(1,Nx-1),-1) + diag(1*sig.*ones(1,Nx-2),-2) + diag(-4*sig.*ones(1,Nx-1),1) + diag(1*sig.*ones(1,Nx-2),2);

M(1,:) = R1;
M(2,:) = R2;
M(Nx-1,:) = RN1;
M(Nx,:) = RN;

M = sparse(M);

count = 10;
tic;
for n = 1:Nt

    % Finite differences.
    ySave(:,2) = inv(M)*ySave(:,1);

    % Bead model.
    X71  = xbSave71(:,:,1);
    Xc71 = xc71(:,n);
    F71 = zeros(3,Nb71);
    U71 = zeros(3,Nb71);

    %% bending forces.
    Fb71 = get_bending_forces(X71);
    F71 = F71 + (Nb71-1).*Fb71;

    %% spring forces.
    Fs71 = get_spring_forces(X71, b071);
    F71  = F71 + calS.*Fs71;

    %% hydrodynamic interactions.
    for p = 1:Nb71
        clear stokeslets71
        stokeslets71 = get_stokeslets(X71, xbSave71(:,p,1),epsilon);
        G71 = reshape(F71',[1,3*Nb71])';
        U71(:,p) = stokeslets71*G71;
    end


    %% update positions and centre of mass.
    xbSave71(:,:,2) = xbSave71(:,:,1) + U71*dt;
    xc71(:,n+1) = mean(xbSave71(:,:,2),2);

    yc(n) = mean(ySave(:,2));

%     if mod(n,viewInt)== 0
%         clf
%         hold on
%         figFD = plot(x,ySave(:,1),'Linewidth',2);
%         figBM = draw_bm_relaxing(X,Xc);
%         figYc = scatter(0.5,yc(1,n),'bx','Linewidth',2);
%         axis([-0.05 1.05 -0.02 0.12])
%         xlabel('x')
%         ylabel('y')
%         box on
%         pause(0.01)
%         hold off
%         disp('fig drawn...')
%     end

    if mod(n,(Nt-1)/10)==0
        fprintf('Script %g perc. complete...',count)
        save('workspace_multiNb_part.mat')
        count = count+10;
    end

    %% replace saved data.
    xbSave71(:,:,1) = xbSave71(:,:,2);

    ySave(:,1) = ySave(:,2);

    %% store points of interest.
    xbMid71(:,n) = xbSave71(:,ceil(Nb71/2),1);    % mid bead
    xbEnd71(:,n) = xbSave71(:,Nb71,1);            % end bead

    yMid(:,n)  = ySave(ceil(Nx/2),1);             % y mid
    yEnd(:,n)  = ySave(Nx,1);                     % y end

end
runtime = toc;

%% Completion.

save('workspace_multiNb.mat')
disp('Workspace saved.')
disp('Script complete.')

%% Function definitions.

function y = intConfig(x)
% Defines the shape of the filament at initial time.
y = 0.6.*(x-0.5).^2;
end

function y = intConfigLegacy(x)
% Defines the shape of the filament at initial time.
% This function was used in the original relaxing rod tests of the bead model.
y = 0.4.*(x-0.5).^2;
end

function fig = draw_bm_relaxing(x,xc)

Nb = size(x,2);

fig1 = scatter3(x(1,:),x(2,:),x(3,:),30,'ro','filled');
axis([-0.1 1.1 -0.1 0.2 0 1])
view(2)
% axis equal
xlabel('x')
ylabel('y')
zlabel('z')
grid on
hold on
lines = zeros(3,100,Nb-1);
for kk = 1:Nb-1
    lines(1,:,kk) = linspace(x(1,kk),x(1,kk+1));
    lines(2,:,kk) = linspace(x(2,kk),x(2,kk+1));
    lines(3,:,kk) = linspace(x(3,kk),x(3,kk+1));
    fig = plot3(squeeze(lines(1,:,kk)),squeeze(lines(2,:,kk)),squeeze(lines(3,:,kk)),'r');
end
scatter3(xc(1),xc(2),xc(3),50, 'rx','Linewidth',2)

end

function F = get_bending_forces(x)
%% GET_BENDING_FORCES
%   Calculates restorative bending forces at Np 3-dimensional beads to keep
%   a chain aligned.
%   INPUTS: x (3xNp matrix of bead positions)
%   OUTPUTS: F (bending forces)

%% Main.

Np = size(x,2);
F = zeros(3,Np);

% bead 1.
bj   = x(:,1) - x(:,2);
bjp1 = x(:,2) - x(:,3);
F(:,1) = bjp1/norm(bj)/norm(bjp1) + dot(bj,bjp1)*(-bj)/norm(bj)^3/norm(bjp1);

% bead 2.
bjn1 = x(:,1) - x(:,2);
bj   = x(:,2) - x(:,3);
bjp1 = x(:,3) - x(:,4);
Fb = (bjn1 - bj)/norm(bjn1)/norm(bj) + dot(bjn1,bj)*(bjn1/norm(bjn1)^3/norm(bj) - bj/norm(bjn1)/norm(bj)^3);
Fc = bjp1/norm(bj)/norm(bjp1) + dot(bj,bjp1)*(-bj)/norm(bj)^3/norm(bjp1);
F(:,2) = Fb+Fc;

% beads 3 to Np-2.
for kk = 3:(Np-2)
    bjn2 = x(:,kk-2) - x(:,kk-1);
    bjn1 = x(:,kk-1) - x(:,kk);
    bj   = x(:,kk)   - x(:,kk+1);
    bjp1 = x(:,kk+1) - x(:,kk+2);
    Fa = -bjn2/norm(bjn2)/norm(bjn1)   + dot(bjn2,bjn1)*bjn1/norm(bjn2)/norm(bjn1)^3;
    Fb = (bjn1-bj)/norm(bjn1)/norm(bj) + dot(bjn1,bj)*(bjn1/norm(bjn1)^3/norm(bj) - bj/norm(bjn1)/norm(bj)^3);
    Fc = bjp1/norm(bj)/norm(bjp1)      + dot(bj,bjp1)*(-bj)/norm(bj)^3/norm(bjp1);
    F(:,kk) = Fa+Fb+Fc;
end

% bead Np-1.
bjn2 = x(:,Np-3) - x(:,Np-2);
bjn1 = x(:,Np-2) - x(:,Np-1);
bj   = x(:,Np-1) - x(:,Np);
Fa = -bjn2/norm(bjn2)/norm(bjn1)   + dot(bjn2,bjn1)*bjn1/norm(bjn2)/norm(bjn1)^3;
Fb = (bjn1-bj)/norm(bjn1)/norm(bj) + dot(bjn1,bj)*(bjn1/norm(bjn1)^3/norm(bj) - bj/norm(bjn1)/norm(bj)^3);
F(:,Np-1) = Fa+Fb;

% bead Np.
bjn2 = x(:,Np-2) - x(:,Np-1);
bjn1 = x(:,Np-1) - x(:,Np);
F(:,Np) = -bjn2/norm(bjn2)/norm(bjn1) + dot(bjn2,bjn1)*bjn1/norm(bjn2)/norm(bjn1)^3;

end

function F = get_spring_forces(x, b0)
%% GET_SPRING_FORCES
%   Calculates the spring forces on a connected chain of beads in 3-dimensional space.
%   INPUTS: X (3xNp matrix of bead locations)
%           b0 (equilibrium distance between beads)
%   OUTPUTS: F (3xNp matrix of forces)

%% Main.

Np = size(x,2);
F  = zeros(3,Np);

for p = 1:Np

    if p > 1
        bpn1  = x(:,p-1) - x(:,p);
        Nbpn1 = norm(bpn1);
        F1    = -2*(Nbpn1 - b0)*(-bpn1/Nbpn1);
    else
        F1 = 0;
    end
    if p < Np
        bp   = x(:,p) - x(:,p+1);
        Nbp  = norm(bp);
        F2   = -2*(Nbp - b0)*(bp/Nbp);
    else
        F2 = 0;
    end

    F(:,p) = F1 + F2;
end

end

function S = get_stokeslets(x, x0, eps)
%GET_REG_STOKESLETS Calculates the regularised stokeslets for a given list of nodes and regularisation parameter.
%   ______________________________________________________________________
%   INPUTS
%   q_nodes: positions of quadrature nodes along the boundary of object.
%   f_nodes: positions of force nodes on the boundary.
%   eps: the regularisation parameter epsilon.
%   mu: the dynamic viscosity of the fluid in question.
%
%   q_nodes and f_nodes must be inputted as a 3xN/3xQ
%   ______________________________________________________________________
%   OUTPUTS
%   S: the block matrix of stokeslets
%   ______________________________________________________________________

eps2 = eps^2;
Q = size(x0,2);

for p = 1:Q

    rx = x(1,:) - x0(1,p);
    ry = x(2,:) - x0(2,p);
    rz = x(3,:) - x0(3,p);

    r2 = rx.^2 + ry.^2 + rz.^2;
    r_eps = (r2 + eps2).^1.5;
    r_num = (r2 + 2*eps2);

    Sxx(p,:) = (r_num + rx.*rx)./r_eps;
    Sxy(p,:) = (rx.*ry)./r_eps;
    Sxz(p,:) = (rx.*rz)./r_eps;

    Syx(p,:) = (ry.*rx)./r_eps;
    Syy(p,:) = (r_num + ry.*ry)./r_eps;
    Syz(p,:) = (ry.*rz)./r_eps;

    Szx(p,:) = (rz.*rx)./r_eps;
    Szy(p,:) = (rz.*ry)./r_eps;
    Szz(p,:) = (r_num + rz.*rz)./r_eps;

end

A1 = horzcat(Sxx, Sxy, Sxz);
A2 = horzcat(Syx, Syy, Syz);
A3 = horzcat(Szx, Szy, Szz);

S = vertcat(A1, A2, A3);

end
