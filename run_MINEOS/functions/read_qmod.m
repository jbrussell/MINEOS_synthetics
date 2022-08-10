% Function to read in the initial model card and plot it to double check
% that everything looks ok
%
% varargin used to take in name of model card for iteration step.
% if no name is given then the model is named init_model
% NJA, 2014

function qmod=read_qmod(QMOD)

isfigure = 0;
% Parameters spcefic to the format of these files
hlines1 = 1;

fid = fopen([QMOD],'r');
A=textscan...
    (fid,'%f %f %f','headerlines',hlines1);

qmod.fname=QMOD;
qmod.z = 6371-A{1};
qmod.rad = A{1}*1000;
qmod.qmu = A{2};
qmod.qkap = A{3};


fclose(fid);

% Choose which card you want to look at

if isfigure
yaxis=[0 350];

figure(1)
clf
plot(qmod.qmu,qmod.z,'-','color','r','linewidth',2);
hold on

plot(qmod.qkap,qmod.z,'-k','linewidth',2);

ylim(yaxis)
xlim([0 1500])
set(gca,'ydir','reverse','fontsize',16)

hleg = legend('Q_{\mu}','Q_{\kappa}');


end
