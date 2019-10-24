
run = 11;
plot3([0 rs_matrix_s(run,5)],[0 rs_matrix_s(run,6)],[0 rs_matrix_s(run,7)],'MarkerSize',5,'Color','b')
hold on
plot3([0 rs_matrix_s(run,8)],[0 rs_matrix_s(run,9)],[0 rs_matrix_s(run,10)],'MarkerSize',5,'Color','r')
axis([-1.2 1.2 -1.2 1.2 -1.2 1.2])
grid('on')
xlabel('x')
ylabel('y')
zlabel('z')

%plot3([0 rs_matrix(1,5)],[0 rs_matrix(1,6)],[0 rs_matrix(1,7)],'MarkerSize',5,'Color','b')
%plot3([0 rs_matrix(1,8)],[0 rs_matrix(1,9)],[0 rs_matrix(1,10)],'MarkerSize',5,'Color','r')