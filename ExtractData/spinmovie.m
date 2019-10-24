%2 - 2.5s 

color1 = [255 81 71]/255;
color2=[63 137 255]/255;
color3 = [242 152 157]/255;
color4=[135 206 250]/255;

i_low = 9900;
i_high = 11000;
i_range = i_low:i_high;

for i=i_range
    plot3([0 s_He(1,i)],[0 s_He(2,i)],[0 s_He(3,i)],'MarkerSize',5,'Color','r')
    hold on
    plot3([0 s_n(1,i)],[0 s_n(2,i)],[0 s_n(3,i)],'MarkerSize',5,'Color','b')
    plot3([1.2 1.2],[0 s_He(2,i)],[0 s_He(3,i)],'MarkerSize',5,'Color',color1)
    plot3([1.2 1.2],[0 s_n(2,i)],[0 s_n(3,i)],'MarkerSize',5,'Color',color2)
    plot3([0 s_He(1,i)],[1.2 1.2],[0 s_He(3,i)],'MarkerSize',5,'Color',color1)
    plot3([0 s_n(1,i)],[1.2 1.2],[0 s_n(3,i)],'MarkerSize',5,'Color',color2)
    plot3([0 s_He(1,i)],[0 s_He(2,i)],[-1.2 -1.2],'MarkerSize',5,'Color',color1)
    plot3([0 s_n(1,i)],[0 s_n(2,i)],[-1.2 -1.2],'MarkerSize',5,'Color',color2)
    plot3(s_He(1,i_low:i),s_He(2,i_low:i),s_He(3,i_low:i),'MarkerSize',5,'Color',color3);
    plot3(s_n(1,i_low:i),s_n(2,i_low:i),s_n(3,i_low:i),'MarkerSize',5,'Color',color4);
    hold off
    axis([-1.2 1.2 -1.2 1.2 -1.2 1.2])
    title(t(1,i),'FontSize',8)
    grid('on')
    xlabel('x')
    ylabel('y')
    zlabel('z')
    pause(0.02)
end

%hold on
%plot3(s_He(1,i_range),s_He(2,i_range),s_He(3,i_range),'MarkerSize',5,'Color','r');
%plot3(s_n(1,i_range),s_n(2,i_range),s_n(3,i_range),'MarkerSize',5,'Color','b');
%hold off
%{
s_He(:,73220)
s_n(:,73220)

He_theta = acos(s_He(1,73220))
n_theta = acos(s_n(1,73220))

He_phi = acos(s_He(2,73220)/sin(He_theta))
n_phi  = acos(s_n(2,73220)/sin(n_theta))
%}