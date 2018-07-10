%% Import workspace and plot figures for relaxing tests.

clear all; clc; clf;
load('workspace_rapid_dynamics.mat')

%% Initial config.
xInt    = x;
yInt    = intConfig(xInt);
ycInt   = yc(1);
yIntLeg = intConfigLeg(xInt);

Nb           = Nb71;
beadInt      = zeros(3,Nb);
beadInt(1,:) = linspace(0,1,Nb);
beadInt(2,:) = intConfig(beadInt(1,:));

beadIntLeg      = zeros(3,Nb);
beadIntLeg(1,:) = linspace(0,1,Nb);
beadIntLeg(2,:) = intConfigLeg(beadIntLeg(1,:));

h1          = subplot(1,3,1);
hold on
box on
intBead     = draw_bm_relaxing(beadInt);
intBeadLeg  = draw_bm_relaxing_leg(beadIntLeg);
lgd1        = legend([intBead,intBeadLeg],'New initial config.','Legacy initial config.','Location','southoutside');

set(lgd1,'Interpreter','latex')
hold off

%% Final config.
h1          = subplot(1,3,2);
beadFin     = xbSave71(:,:,2);

hold on
box on
finBead     = draw_bm_relaxing(beadFin);
% finBeadLeg  = draw_bm_relaxing_leg(beadIntLeg);

hold off


%% Fast dynamics insert

h1 = subplot(1,3,3);
hold on
box on
BMc   = plot(t,xc71(2,2:end),'r','LineWidth',1);
% xlabel('time')
% axis([0 1.4e-5 0.03423 0.03432])
% ylabel('$y$','Interpreter','latex')
hold off


%% Resize figure window and save pdf
% set(h1,'pos',[100 100 520 550])
% set(h2,'pos',[900 100 860 550])
% pos = get(gcf,'Position');
% set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3),pos(4)])
% print(gcf,'compPlots71Mat','-dpdf','-r0')

%% Functions

function fig = draw_bm_relaxing(x)
% edited to remove centre of mass plot.

Nb = size(x,2);

fig1 = scatter3(x(1,:),x(2,:),x(3,:),30,'bo','filled');
% axis([-0.1 1.1 -0.1 0.2 0 1])
view(2)
% axis equal
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
zlabel('$z$','Interpreter','latex')
grid on
hold on
lines = zeros(3,100,Nb-1);
for kk = 1:Nb-1
    lines(1,:,kk) = linspace(x(1,kk),x(1,kk+1));
    lines(2,:,kk) = linspace(x(2,kk),x(2,kk+1));
    lines(3,:,kk) = linspace(x(3,kk),x(3,kk+1));
    fig = plot3(squeeze(lines(1,:,kk)),squeeze(lines(2,:,kk)),squeeze(lines(3,:,kk)),'b');
end
% scatter3(xc(1),xc(2),xc(3),50, 'kx','Linewidth',2)

end

function fig = draw_bm_relaxing_leg(x)
% edited to remove centre of mass plot.
% plots beads a different colour compared to default plotting function.

Nb = size(x,2);

fig1 = scatter3(x(1,:),x(2,:),x(3,:),30,'ro','filled');
% axis([-0.1 1.1 -0.1 0.2 0 1])
view(2)
% axis equal
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
zlabel('$z$','Interpreter','latex')
grid on
hold on
lines = zeros(3,100,Nb-1);
for kk = 1:Nb-1
    lines(1,:,kk) = linspace(x(1,kk),x(1,kk+1));
    lines(2,:,kk) = linspace(x(2,kk),x(2,kk+1));
    lines(3,:,kk) = linspace(x(3,kk),x(3,kk+1));
    fig = plot3(squeeze(lines(1,:,kk)),squeeze(lines(2,:,kk)),squeeze(lines(3,:,kk)),'r');
end
% scatter3(xc(1),xc(2),xc(3),50, 'kx','Linewidth',2)

end

function y = intConfig(x)
% Defines the shape of the filament at initial time.
y = 0.6.*(x-0.5).^2;
end

function y = intConfigLeg(x)
% Defines the shape of the filament at initial time.
% This function was used in the original relaxing rod tests of the bead model.
y = 0.4.*(x-0.5).^2;
end
