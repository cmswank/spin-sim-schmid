

%% Parameters import %%
run_num = 3434;
simType = 0;
runTime = 0.009;
Bin = 101;
Event = 1;
fdir = '/home/ezra/spin-sim-SD/ExtractData/';
datafile = strcat(fdir,'OutputHeatflush',num2str(run_num),'.dat');


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
t = squeeze(tlarge(1,:));


%Constants etc.
k = 1.38064852e-23;         	% Boltzmann const.
gamma_n = -1.8324e8;       		% neutron gyromagnetic ratio, rad s^-1 T^-1
gamma_3 = -2.037947093e8;   	% He3 gyromagnetic ratio, rad s^-1 T^-1
M_3 = 2.2*3.016*1.66054e-27;  	% effective mass of He3 in super fluid
M_n = 1.008664*1.66054e-27;		% neutron mass
L_x = 0.076;                	% trap x dim /m (H_0 along x by convention)
L_y = 0.102;               		% trap y dim /m
L_z = 0.400;                  	% trap z dim /m

% n and He3 ez import
if simType == 1
	s_He = s;
end
if simType == 0
	s_n = s;
end

s_d = s_n - s_He;

cros = cross(s_He,s_n,1);
pangle = asin(cros(1,:,:)./(sqrt((1-s_n(1,:,:).^2).*(1-s_He(1,:,:).^2)))); % angle of projection on y z plane between He and n
angle = acos(dot(s_n,s_He,1)); % total angle between He and n 
signal = 1-dot(s_n,s_He,1);
plot(t,pangle(1,:));
%plot(t,anglerun50(1,:,1));
hold on
plot(t,angle(1,:,1));
legend('pangle','angle')
hold on
plot([0 runTime], [0.8 0.8])
%plot([0 runTime], [-0.8 -0.8])
xlabel('t /s')
ylabel('angle /rad')
%figure;
%plot(t, signal(1,:,1));
%xlabel('t /s')
%ylabel('1-cos(angle)')
%plot(t,pangle(1,:,1));
%xlabel('t /s')
%ylabel('phi /rad')
%figure;
%plot(t,sqrt(s_n(2,:).^2+s_n(3,:).^2))
