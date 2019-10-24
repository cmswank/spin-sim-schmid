
frame = 73220;
s_He(:,frame);
s_n(:,frame);

He_theta_50 = acos(s_He(1,frame));
n_theta_50 = acos(s_n(1,frame));

He_phi_50 = acos(s_He(2,frame)/sin(He_theta_50));
n_phi_50  = acos(s_n(2,frame)/sin(n_theta_50));
angle_50 = acos(dot(s_n(:,frame),s_He(:,frame),1));

He_theta = 1.7939;
n_theta = 1.6686;

He_phi = 0.5176;
n_phi  = 0.4710;

s_n_i = [cos(n_theta) sin(n_theta)*cos(n_phi) sin(n_theta)*sin(n_phi)]
s_He_i = [cos(He_theta) sin(He_theta)*cos(He_phi) -sin(He_theta)*sin(He_phi)]
angle_i = acos(dot(s_n_i,s_He_i,2))