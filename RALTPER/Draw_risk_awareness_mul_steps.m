function Draw_risk_awareness_mul_steps(xs,T,N,G,e_g,X_set,U_set,A_set,F_set,b_set,fig)
mpc_N = size(X_set,3);

Rx = @(x)[cos(x),sin(x);-sin(x),cos(x)];

set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize', 12)

line_width = 1.5;
fontsize_labels = 14;


figure(fig)

set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');
set(gcf,'Units','normalized','OuterPosition',[0 0 0.55 1]);

for k = 1 : mpc_N
    % hold on
    axis([-5 30 -1 4])
    fill([0,0,30,30],[0,-4,-4,0],'k');
    hold on
    fill([0,30,30,0],[3,3,4,4],'k');
    h_t = 0.5; w_t=0.3; % triangle parameters
    x1 = xs(1); y1 = xs(2); th1 = xs(3);
    x1_tri = [ x1+h_t*cos(th1), x1+(w_t/2)*cos((pi/2)-th1), x1-(w_t/2)*cos((pi/2)-th1)];%,x1+(h_t/3)*cos(th1)];
    y1_tri = [ y1+h_t*sin(th1), y1-(w_t/2)*sin((pi/2)-th1), y1+(w_t/2)*sin((pi/2)-th1)];%,y1+(h_t/3)*sin(th1)];
    fill(x1_tri, y1_tri, 'g'); % plot reference state
    X0 = X_set(:,:,k);
    for j = 1:1
        tmp_off = X0(j,1:2)';
        tmp_alpha = X0(j,3);
        tmp_R = Rx(tmp_alpha);
        hyperplane_plot(G*tmp_R,e_g+G*tmp_R*tmp_off,0.2);
    end
    
    t = linspace(0, 2*pi, 1000);
    for j = 1:1
        F = F_set(:,:,j,k);
        A = A_set(:,:,j,k);
        b = b_set(:,j,k);
        % 计算椭圆上的点
        ellipse_points = [1/sqrt(F(1,1)) * cos(t); 1/sqrt(F(2,2)) * sin(t)];
        rotated_points = A' * ellipse_points;
        X = rotated_points(1, :) + b(1);
        Y = rotated_points(2, :) + b(2);
        fill(X,Y,'r','FaceAlpha',0.5);
    end
    
    plot(X0(:,1),X0(:,2),'-*r','linewidth',line_width)
    
    hold off
    ylabel('$y$-position (m)','interpreter','latex','FontSize',fontsize_labels)
    xlabel('$x$-position (m)','interpreter','latex','FontSize',fontsize_labels)
    pause(0.1)
    box on;
    grid on
    drawnow
    FF(k) = getframe(gcf); % to get the current frame
end
close(gcf)

video = VideoWriter('exp.avi','Motion JPEG AVI');
video.FrameRate = 5; 
open(video)
writeVideo(video,FF)
close (video)
figure()
axis([-5 30 -1 4])
fill([0,0,30,30],[0,-4,-4,0],'k');
hold on
fill([0,30,30,0],[3,3,4,4],'k');
for j = 1:mpc_N
    tmp_off = X_set(1,1:2,j)';
    tmp_alpha = X_set(1,3,j);
    tmp_R = Rx(tmp_alpha);
    hyperplane_plot(G*tmp_R,e_g+G*tmp_R*tmp_off,0.2);
end


figure()
subplot(211)
plot((1:mpc_N)*T,U_set(:,1));
ylabel('v (rad/s)')
grid on
subplot(212)
plot((1:mpc_N)*T,U_set(:,2));
xlabel('time (seconds)')
ylabel('\omega (rad/s)')
grid on

end