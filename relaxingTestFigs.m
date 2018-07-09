%% Import workspace and plot figures for relaxing tests.

clear all; clc; clf;
load('workspace_multiNb.mat')

%% Initial config.
xInt  = x;
yInt  = intConfig(xInt);
ycInt = yc(1); 

beadInt11      = zeros(3,Nb11);
beadInt11(1,:) = linspace(0,1,Nb11);
beadInt11(2,:) = intConfig(beadInt11(1,:));
beadcInt11     = xc11(:,1);

beadInt71      = zeros(3,Nb71);
beadInt71(1,:) = linspace(0,1,Nb71);
beadInt71(2,:) = intConfig(beadInt71(1,:));
beadcInt71     = xc71(:,1);

h1 = figure(1);
subplot(2,2,1)
hold on 
box on 
FDInt1  = plot(xInt,yInt,'b', 'Linewidth',2);
FDcInt1 = scatter(0.5,ycInt,'bo','filled');
BMInt1  = draw_bm_relaxing(beadInt11,beadcInt11);
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
title('Initial configuration')
axis([0 1 -0.01 0.11])
legend('off')
hold off

h2 = figure(2);
subplot(2,3,1)
hold on 
box on 
FDInt2  = plot(xInt,yInt, 'Linewidth',2);
FDcInt2 = scatter(0.5,ycInt,'bo','filled');
BMInt2  = draw_bm_relaxing(beadInt71,beadcInt71);
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
title('Initial configuration')
axis([0 1 -0.01 0.11])
legend('off')
hold off


%% Final config.
xFin  = x;
yFin  = ySave(:,2);
ycFin = yc(Nt);

beadFin11  = xbSave11(:,:,2);
beadcFin11 = xc11(:,Nt);
beadFin71  = xbSave71(:,:,2);
beadcFin71 = xc71(:,Nt);

h1 = figure(1);
subplot(2,2,2)
hold on
box on
FDFin1  = plot(xFin,yFin,'b','LineWidth',2);
FDcInt1 = scatter(0.5,ycFin,'bo','filled');
BMFin1  = draw_bm_relaxing(beadFin11,beadcFin11);
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
title('Final configuration')
axis([0 1 -0.01 0.11])
legend('off')
hold off

h2 = figure(2);
subplot(2,3,2)
hold on
box on
FDFin2  = plot(xFin,yFin,'LineWidth',2);
FDcInt2 = scatter(0.5,ycFin,'bo','filled');
BMFin2  = draw_bm_relaxing(beadFin71,beadcFin71);
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
title('Final configuration')
axis([0 1 -0.01 0.11])
legend('off')
hold off

%% Mid/End point travel.

h1 = figure(1);
subplot(2,2,3)
hold on
box on
BMc11   = plot([1:Nt],xc11(2,2:end),'--r','LineWidth',1);
BMMid11 = plot([1:Nt],xbMid11(2,:),'r');
BMEnd11 = plot([1:Nt],xbEnd11(2,:),'r');
FDc   = plot([1:Nt],yc,'--b','LineWidth',1);
FDMid = plot([1:Nt],yMid,'b');
FDEnd = plot([1:Nt],yEnd,'b');
xlabel('time')
ylabel('$y$','Interpreter','latex')
axis tight
title('Mid/end point travel')
lgd11 = legend([BMc11,FDc],{'BM centre of mass','EHD model centre of mass'},'Location','southoutside');
set(gca,'xticklabel',[])
hold off

h2 = figure(2);
subplot(2,3,4)
hold on
box on
BMc71   = plot([1:Nt],xc71(2,2:end),'--r','LineWidth',1);
BMMid71 = plot([1:Nt],xbMid71(2,:),'r');
BMEnd71 = plot([1:Nt],xbEnd71(2,:),'r');
FDc   = plot([1:Nt],yc,'--b','LineWidth',1);
FDMid = plot([1:Nt],yMid,'b');
FDEnd = plot([1:Nt],yEnd,'b');
xlabel('time')
ylabel('$y$','Interpreter','latex')
axis tight
title('Mid/end point travel')
lgd12 = legend([BMc71,FDc],{'BM centre of mass','EHD model centre of mass'},'Location','southoutside');
set(gca,'xticklabel',[])
hold off


%% CoM travel.

h1 = figure(1);
subplot(2,2,4)
hold on
box on
BMc11 = plot([1:Nt],xc11(2,2:end),'r','LineWidth',1);
FDc   = plot([1:Nt],yc,'b','LineWidth',1);
xlabel('time')
ylabel('$y$','Interpreter','latex')
axis([0 5e5 0.033 0.041])
title('Centre of mass travel')
set(gca,'xticklabel',[])
lgd21 = legend('Bead model','EHD model','Location','southoutside');
hold off

h2 = figure(2);
subplot(2,3,5)
hold on
box on
BMc71 = plot([1:Nt],xc71(2,2:end),'r','LineWidth',1);
FDc   = plot([1:Nt],yc,'b','LineWidth',1);
xlabel('time')
ylabel('$y$','Interpreter','latex')
axis([0 5e5 0.033 0.041])
title('Centre of mass travel')
set(gca,'xticklabel',[])
lgd22 = legend('Bead model','EHD model','Location','southoutside');
hold off

%% Fast dynamics insert
h2 = figure(2);
subplot(2,3,6)
hold on
box on
BMc   = plot(t,xc71(2,2:end),'r','LineWidth',1);
% xlabel('time')
axis([0 1.4e-5 0.03423 0.03432])
% ylabel('$y$','Interpreter','latex')
hold off


%% Resize figure window and save pdf
set(h1,'pos',[100 100 520 550])
set(h2,'pos',[900 100 860 550])
% pos = get(gcf,'Position');
% set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3),pos(4)])
% print(gcf,'compPlots71Mat','-dpdf','-r0')

%% Functions 

function fig = draw_bm_relaxing(x,xc)

Nb = size(x,2);

fig1 = scatter3(x(1,:),x(2,:),x(3,:),30,'ro','filled');
% axis([-0.1 1.1 -0.1 0.2 0 1]) 
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
scatter3(xc(1),xc(2),xc(3),50, 'kx','Linewidth',2)

end

function y = intConfig(x)
% Defines the shape of the filament at initial time.
y = 0.4.*(x-0.5).^2;
end
