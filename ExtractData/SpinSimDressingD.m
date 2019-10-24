
%% Parameters import %%
run_num = 702;
double_import = 1;

fdir = '/home/ezra/spin-sim-SD/ExtractData/';
%fdir = 'E:\Caletch Spin Dressing\SD-runs\';

pfn =  strcat(fdir,'run',num2str(run_num),'_pars'); 	%parameter file name
datafile = strcat(fdir,'OutputHeatflush',num2str(run_num),'.dat');

pfile = fopen(pfn);
pars = textscan(pfile, '%[^= ]%*[= ]%s', 'CommentStyle', '//', 'whitespace','\n', 'delimiter', ';');
fclose(pfile);

simType = str2double(pars{2}(strcmp(pars{1}, 'simType')));
runTime = str2double(pars{2}(strcmp(pars{1}, 'runTime')));
nBins = str2double(pars{2}(strcmp(pars{1}, 'nBins')));
nEntries = str2double(pars{2}(strcmp(pars{1}, 'nEntries')));

T = str2double(pars{2}(strcmp(pars{1}, 'temperature')));
diff_coeff = str2double(pars{2}(strcmp(pars{1}, 'diffusion')));
gravity = str2double(pars{2}(strcmp(pars{1}, 'gravity')));
spinSet_s = pars{2}(strcmp(pars{1}, 'spinSet'));
spinSet = str2double(strsplit(spinSet_s{1},','));

B0_s = pars{2}(strcmp(pars{1}, 'B0'));
B0 = str2double(strsplit(B0_s{1},','));
E0_s = pars{2}(strcmp(pars{1}, 'E0'));
E0 = str2double(strsplit(E0_s{1},','));
%gradB_s = pars{2}(strcmp(pars{1}, 'UniformG'));
%gradB = str2double(strsplit(gradB_s{1},','));

SD_s = pars{2}(strcmp(pars{1}, 'SpinDressing'));
SD = str2double(strsplit(SD_s{1},','));
w_rf = SD(2);
phi = SD(3);
if SD(1) == 5
	fm = SD(4);
	dw1 = SD(5);
	dw2 = SD(6);
	deltat1 = SD(7);
    deltat2 = SD(8);
end
B_rf = str2double(strrep( pars{2}(strcmp(pars{1}, 'BzAdd')),'"',''));

Bin = nBins + 1;
Event = nEntries;

%% import data %%
fileID = fopen(datafile);
A = fread(fileID, 'double');
fclose(fileID);

B = reshape(A, 10, Bin, Event);
% matrices (time (Bin), particles (Event))
s = squeeze(B(1:3,:,:));
sx = squeeze(B(1,:,:));
sy = squeeze(B(2,:,:));
sz = squeeze(B(3,:,:));
x =  squeeze(B(4,:,:));
y =  squeeze(B(5,:,:));
z =  squeeze(B(6,:,:));
vx = squeeze(B(7,:,:));
vy = squeeze(B(8,:,:));
vz = squeeze(B(9,:,:));
tlarge = squeeze(B(10,:,:)); 
t = squeeze(tlarge(:,1));


%Constants etc.
k = 1.38064852e-23;         	% Boltzmann const.
gamma_n = -1.83247185e8;       	% neutron gyromagnetic ratio, rad s^-1 T^-1, extracted from sim code
gamma_3 = -2.037894730e8;   	% He3 gyromagnetic ratio, rad s^-1 T^-1
M_3 = 2.2*3.016*1.66054e-27;  	% effective mass of He3 in super fluid
M_n = 1.008664*1.66054e-27;		% neutron mass
L_x = 0.076;                	% trap x dim /m (H_0 along x by convention)
L_y = 0.102;               		% trap y dim /m
L_z = 0.400;                  	% trap z dim /m
D = 1.6e-4/T^7;             	% mass diffusion coefficient of He3 m^2 s^-1
v_xsqrd_3 = k*T/M_3;            % average value of (x comp of vel squared) of He3
v_rms_3 = sqrt(3*v_xsqrd_3);    % rms speed of He3
tau_c_3 = D/sqrt(v_xsqrd_3);    % mean free period of He3

% n and He3 ez import
if simType == 1
	s_He = s;
else
	s_n = s;
end

if double_import == 1
    datafile = strcat(fdir,'OutputHeatflush',num2str(run_num+1),'.dat');
    fileID = fopen(datafile);
    A = fread(fileID, 'double');
    fclose(fileID);
    B = reshape(A, 10, Bin, Event);
    s = squeeze(B(1:3,:,:));
    if simType == 1
        s_n = s;
    else
        s_He = s;
    end
end

s_d = s_n - s_He;


%shoun't be n and He below but are because I'm lazy
s_n = sum(s_n,3)./nEntries;
s_He = sum(s_He,3)./nEntries;
cros = cross(s_n,s_He,1);
pangle = asin(cros(1,:)./(sqrt((1-s_n(1,:,1).^2).*(1-s_He(1,:,1).^2))));
fc = 100; % Cut off frequency
fs = nBins/runTime; % Sampling rate
[b,a] = butter(6,fc/(fs/2)); % Butterworth filter of order 6
%angle_f = filter(b,a,angle(1,:,1)); % Will be the filtered signal
pangle_f = filter(b,a,pangle(1,:));
plot(t,pangle_f);
hold on;
plot(t,conv(pangle, ones(501,1)/501, 'same'));
xlabel('t /s')
ylabel('angle /rad')
legend('angle','angle (moving average)')
m = t\pangle.';
plot(t,t*m);



%{
cros = cross(s_He,s_n,1);
pangle = asin(cros(1,:,:)./(sqrt((1-s_n(1,:,:).^2).*(1-s_He(1,:,:).^2)))); % angle of projection on y z plane between He and n
angle = acos(dot(s_n,s_He,1)); % total angle between He and n 
signal = 1-dot(s_n,s_He,1);

fc = 100; % Cut off frequency
fs = nBins/runTime; % Sampling rate
[b,a] = butter(6,fc/(fs/2)); % Butterworth filter of order 6
angle_f = filter(b,a,angle(1,:,1)); % Will be the filtered signal
pangle_f = filter(b,a,pangle(1,:));
%plot(t,pangle(1,:,1));
%hold on
%plot(t,angle(1,:,1));
plot(t,pangle_f);
hold on
plot(t,angle_f);
%legend('pangle','angle','pangle_f','angle_f')
legend('angle','projected angle')
plot([0 runTime], [0.8 0.8])
plot(t,ModulationFunc(2*pi*t+pi/2,51,0.726))
xlabel('t /s')
ylabel('angle /rad')
%plot([0 runTime], [-0.8 -0.8])
%}




%{
plot(t,pangle(1,:,1));
hold on
plot(t,angle(1,:,1));
legend('pangle','angle')
plot([0 runTime], [0.8 0.8])
%plot([0 runTime], [-0.8 -0.8])
xlabel('t /s')
ylabel('\Delta\theta /rad')
figure;

%plot(t, signal(1,:,1));
%xlabel('t /s')
%ylabel('1-cos(angle)')
%plot(t,pangle(1,:,1));
%xlabel('t /s')
%ylabel('phi /rad')
%figure;
%plot(t,sqrt(s_n(2,:).^2+s_n(3,:).^2))

pangle_n = atan(s_n(2,:)./s_n(3,:));
pangle_He = atan(s_He(2,:)./s_He(3,:));
plot(t, pangle_n)
hold on
plot(t, pangle_He)
legend('n','He')
xlabel('t /s')
ylabel('\phi /rad')
figure;
%}


%{
ax1 = subplot(3,1,1);
ax2 = subplot(3,1,2);
ax3 = subplot(3,1,3);
plot(ax1,t(:),s_n(1,:),'b');
hold on
plot(ax1,t(:),s_He(1,:),'r');
plot(ax2,t(:),s_n(2,:),'b');
plot(ax2,t(:),s_He(2,:),'r');
plot(ax3,t(:),s_n(3,:),'b');
plot(ax3,t(:),s_He(3,:),'r');
xlabel(ax1,'t /s')
xlabel(ax2,'t /s')
xlabel(ax3,'t /s')
figure;

ax1 = subplot(3,1,1);
ax2 = subplot(3,1,2);
ax3 = subplot(3,1,3);
plot(ax1,t(:),s_d(1,:),'r');
hold on
plot(ax2,t(:),s_d(2,:),'r');
plot(ax3,t(:),s_d(3,:),'r');
xlabel(ax1,'t /s')
xlabel(ax2,'t /s')
xlabel(ax3,'t /s')


ax1 = subplot(3,1,1); % top subplot
ax2 = subplot(3,1,2); % bottom subplot
ax3 = subplot(3,1,3);
plot(ax1,t,s_n(2,:));
hold on
plot(ax2,t,s_He(2,:),'r')
plot(ax3,t,s_d(2,:),'g');
title(ax1,'neutrons')
title(ax2,'He3')
title(ax3,'Difference')
xlabel(ax1,'t /s')
xlabel(ax2,'t /s')
xlabel(ax3,'t /s')
ylabel(ax1,'\Sigma S_{y}')
ylabel(ax2,'\Sigma S_{y}')
ylabel(ax3,'\Sigma S_{y,n}-\Sigma S_{y,He_3}')
%}
