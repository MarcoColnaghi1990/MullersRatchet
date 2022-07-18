clf;hold on
ax1=subplot(1,1,1);

s = 10.^(-4:.1:-1);
U = 10.^(-4:.1:-1);

N = 10^5;

n_0 = N *exp(-U./s');

surf(s,U,n_0);

set(gca,'xscale','log','yscale','log')
view([-27.5 30.8]);
colorbar;
% ylim([10^-5, 0.1])
% xlim([10^-5, 0.1])

ax1.Colorbar.Position(1) = 0.9;
ax1.Position(3)=0.7;
ax1.Position(2)=0.15;
set(gca,'FontName','Lucida Bright','Fontsize',12)
ylabel('Strength of selection (\its\rm)','Rotation',-31.5);
ax1.YLabel.Position(2)=0.005;
ax1.YLabel.Position(3)=-20000;
xlabel('Genome-wide mutation rate (\itU\rm)','Rotation',8);
zlabel('Equilibrium LLC size (\itn\rm_0)');

zticks([0,20000,40000,60000,80000,100000])
zticklabels({'0','20000','40000','60000','80000','100000'})