

pfn = '/data1/cmswank/spin-sim-xliu/run30_pars'; 	%parameter file name
datafile = '/data1/cmswank/spin-sim-xliu/ExtractData/SpinDressingCrossTerm_71.dat';

%{
params = fopen(pfn);
for i = 0:7
    tline = fgetl(params);
end
Bin_line = fgetl(params);
Bin = int16(str2double(Bin_line(9:length(Bin_line)-1)) +1);
Event_line = fgetl(params);
Event = int16(str2double(Event_line(12:length(Event_line)-1)));
fclose(params);

%datafile = 'E:\Caletch Spin Dressing\OutputHeatflush20.dat';
%}

Bin = 301;
Event = 1000;

% import data
fileID = fopen(datafile);
A = fread(fileID, 'double');
fclose(fileID);

B = reshape(A, 10, Bin, Event);
% the following are matrices (time (Bin+1), particles (Event-1))
sx = squeeze(B(1,:,2:Event));
sy = squeeze(B(2,:,2:Event));
sz = squeeze(B(3,:,2:Event));
x =  squeeze(B(4,:,2:Event));
y =  squeeze(B(5,:,2:Event));
z = squeeze(B(6,:,2:Event));
vx =  squeeze(B(7,:,2:Event));
vy =  squeeze(B(8,:,2:Event));
vz = squeeze(B(9,:,2:Event));
tlarge = squeeze(B(10,:,2:Event)); 
t = squeeze(tlarge(:,1));

% 
% %McGregor
% k = 1.38064852e-23;         % Boltzmann const.
% gamma_n = -1.8324e8;        % neutron gyromagnetic ratio, rad s^-1 T^-1
% gamma_3 = -2.037947093e8;   % He3 gyromagnetic ratio, rad s^-1 T^-1
% M = 2.2*3.016*1.66054e-27;  % effective mass of He3 in super fluid 
% B0 = 3e-6;                  % mag flux density in x dir, T
% gradB = [ 0 0 5e-8 ; 0 0 0 ; 0 0 -1e-6 ]; % grad of B, gradB(i,j) = derivative in j direction of i comp of B field
% L_x = 0.076;                % trap x dim /m (H_0 along x by convention)
% L_y = 0.102;                % trap y dim /m
% L_z = 0.4;                  % trap z dim /m
% T = 0.450;                  % temp /K
% D = 1.6e-4/T^7;             % mass diffusion coefficient of He3 m^2 s^-1
% v_xsqrd = k*T/M;            % average value of (x comp of vel squared)
% v_rms = sqrt(3*v_xsqrd);    % rms speed of He3
% tau_c = D/sqrt(v_xsqrd);    % mean free period
% 
% gradSqr = sum(abs(gradB).^2,2); % vector of (mod(grad(B_i)))^2
% T_1 = ( (gradSqr(2) + gradSqr(3))*v_xsqrd*tau_c/(B0^2*(1+(B0*gamma_3)^2*tau_c^2)))^-1; % T_1 relaxation time
% T_2 = (0.5*T_1^-1 + gamma_3^2*L_z^4*gradB(1,3)^2/(120*D))^-1;  % T_2 relaxation time
% 
% %22 up, 21 down
% 
% %s_up = squeeze(B(1:3,:,2:Event));
% %s_down = squeeze(B(1:3,:,2:Event));
% 
% %ucd = cross(s_down,s_up,1);
% %dwd = asin(ucd(1,:,:)./(sqrt((1-s_down(1,:,:).^2).*(1-s_up(1,:,:).^2))));
% %{
% mean(dwd(1,Bin,:),3)
% plot(t,mean(dwd(1,:,:),3))
% hold on
% %plot(t,t*(1.044e-3))
% legend('simulation','theory')
% title('Parallel antiparallel phase difference')
% xlabel('t /s')
% ylabel('\Delta \phi')
% figure;
% %}
% %data analysis
% 
% %SY_n = mean(sy,2);  %10
% %SY_3 = mean(sy,2); %11
% SY_d = SY_n - SY_3;
% 
% s_n = squeeze(B(1:3,:,2:Event));
% %s_3 = squeeze(B(1:3,:,2:Event));
% 
% ucd = cross(s_n,s_3,1); % up cross down
% dwd = asin(ucd(1,:,:)./(sqrt((1-s_n(1,:,:).^2).*(1-s_3(1,:,:).^2))));
% 
% plot(t,mean(dwd(1,:,:),3));
% figure;
% 
% nd3 = acos(dot(s_n,s_3,1));
% plot(t,nd3(1,:,1));
% figure;
% 
% 
% ax1 = subplot(3,1,1); % top subplot
% ax2 = subplot(3,1,2); % bottom subplot
% ax3 = subplot(3,1,3);
% plot(ax1,t,SY_n);
% hold on
% plot(ax2,t,SY_3,'r')
% plot(ax3,t,SY_d,'g');
% title(ax1,'neutrons')
% title(ax2,'He_3')
% title(ax3,'Difference')
% xlabel(ax1,'t /s')
% xlabel(ax2,'t /s')
% xlabel(ax3,'t /s')
% ylabel(ax1,'\Sigma S_{y}')
% ylabel(ax2,'\Sigma S_{y}')
% ylabel(ax3,'\Sigma S_{y,n}-\Sigma S_{y,He_3}')
% 
% 
% %{
% 
% % T_1 and T_2
% SX = mean(sx,2);
% SY = mean(sy,2);
% SZ = mean(sz,2);
% 
% SX_i = SX(1)*exp(-t/T_1);
% SY_i = SY(1)*exp(-t/T_2);
% SZ_i = SZ(1)*exp(-t/T_2);
% 
% plot(t,SX)
% hold on
% plot(t,SX_i)
% legend('simulation','theory')
% title('T_1 Relaxation')
% xlabel('t /s')
% ylabel('\Sigma S_{x}')
% hold off
% figure;
% 
% 
% plot(t,SY)
% hold on
% plot(t,SY_i)
% legend('simulation','theory')
% title('T_2 Relaxation')
% xlabel('t /s')
% ylabel('\Sigma S_{y}')
% %plot(t,SZ)
% %plot(t,SZ_i)
% legend('show')
% hold off
% 
% %}
% 

