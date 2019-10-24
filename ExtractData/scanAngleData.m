fdir = '/home/ezra/spin-sim-SD/ExtractData/';


frn = 23000;     % first run number
%{
alpha_steps = 16;
beta_steps = 20;
gamma_steps = 8;
alpha_max = 9*pi/50.0;
alpha_min = -8*pi/50.0;
beta_min = 6*pi/100.0;
beta_max = 0*pi/100.0;
gamma_max = -39*pi/100.0;
gamma_min = -41*pi/100.0;

sep = 0.8;

Bin = 201;
Event = 1;
%}

frn = 1000;     % first run number

alpha_steps = 20;
beta_steps = 20;
gamma_steps = 20;
alpha_max = pi/2;
alpha_min = -pi/2;
beta_max = pi/2;
gamma_max = -pi/2;
gamma_min = -pi/2;

sep = 0.8;

Bin = 101;
Event = 1;

s_n_i = [0; cos(0.5*sep); -sin(0.5*sep)];	% initial neutron spin
s_He_i = [0; cos(0.5*sep); sin(0.5*sep)];	% initial He spin

run = frn;
range_matrix = zeros(alpha_steps*beta_steps*gamma_steps, 2);
prange_matrix = zeros(alpha_steps*beta_steps*gamma_steps, 2);
spin_matrix = zeros(alpha_steps*beta_steps*gamma_steps, 7);
rs_matrix = zeros(alpha_steps*beta_steps*gamma_steps, 13);
prs_matrix = zeros(alpha_steps*beta_steps*gamma_steps, 13);

for i = 0:(alpha_steps-1)
    alpha = (alpha_max - alpha_min)*i/alpha_steps + alpha_min;
	Rx = [[1,0,0]; [0, cos(alpha), -sin(alpha)]; [0, sin(alpha), cos(alpha)]];
    disp(i)
    for j = 0:(beta_steps-1)
        beta = (beta_max-beta_min)*j/beta_steps + beta_min;
        Ry = [[cos(beta), 0, sin(beta)];[0, 1 , 0]; [-sin(beta), 0, cos(beta)]];
        for k = 0:(gamma_steps-1)
            gamma = (gamma_max-gamma_min)*k/gamma_steps + gamma_min;
            Rz = [[cos(gamma), -sin(gamma), 0]; [sin(gamma), cos(gamma), 0]; [0, 0,  1]];
            
            s_n_r = Ry*s_n_i;
			s_He_r = Ry*s_He_i;
			s_n_r = Rx*s_n_r;
			s_He_r = Rx*s_He_r;
			s_n_r = Rz*s_n_r;
			s_He_r = Rz*s_He_r;
    
            datafile = strcat(fdir,'OutputHeatflush',num2str(run),'.dat');
            fileID = fopen(datafile);
            A = fread(fileID, 'double');
            fclose(fileID);
            B = reshape(A, 10, Bin, Event);
            s_n = squeeze(B(1:3,:,:));
            datafile = strcat(fdir,'OutputHeatflush',num2str(run+1),'.dat');
            fileID = fopen(datafile);
            A = fread(fileID, 'double');
            fclose(fileID);
            B = reshape(A, 10, Bin, Event);
            s_He = squeeze(B(1:3,:,:));
            angle = acos(dot(s_n,s_He,1));
            cros = cross(s_He,s_n,1);
            pangle = abs(asin(cros(1,:,:)./(sqrt((1-s_n(1,:,:).^2).*(1-s_He(1,:,:).^2)))));
            min_ang = min(angle(1,:,1));
            max_ang = max(angle(1,:,1));
            ang_range = max_ang - min_ang;
            min_pang = min(angle(1,:,1));
            max_pang = max(pangle(1,:,1));
            pang_range = max_pang - min_pang;
            range_matrix(1+(run-frn)/2, :) = [run ang_range];
            prange_matrix(1+(run-frn)/2, :) = [run pang_range];
            spin_matrix(1+(run-frn)/2, 2:7) = [s_n_r.' s_He_r.'];
            rs_matrix(1+(run-frn)/2, :) = [run ang_range min_ang max_ang s_n_r.' s_He_r.' i j k];
            prs_matrix(1+(run-frn)/2, :) = [run pang_range min_pang max_pang s_n_r.' s_He_r.' i j k];
            run = run + 2;
        end
    end
end

%fileID = fopen('spinAngularRangeData.txt','w');
%fprintf(fileID, '%s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\n', 'run#','ang_range','min_ang','max_ang','s_n_x','s_n_y','s_n_z','s_He_x','s_He_y','s_He_z','x_step','y_step','z_step');
%fprintf(fileID,'%i\t %1.6e\t %1.6e\t %1.6e\t %1.6e\t %1.6e\t %1.6e\t %1.6e\t %1.6e\t %1.6e\t %i\t %i\t %i\n',rs_matrix(:,:));
%fclose(fileID);

[Y,sortlist] = sort(rs_matrix(:,2));
rs_matrix_s = rs_matrix(sortlist,:);

rs_matrix_s(1:10,:);

[Y,sortlist] = sort(prs_matrix(:,2));
prs_matrix_s = prs_matrix(sortlist,:);
