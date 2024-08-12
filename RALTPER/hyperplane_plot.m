function hyperplane_plot(G,g,color,FaceAlpha,flag)

if nargin < 2
    error('函数hyperplane_plot：输入参数过少');
elseif nargin == 2
    color = 'b';
    flag = 0;
    FaceAlpha = 0.3;
elseif nargin == 3
    flag = 0;
    FaceAlpha = 0.3;
elseif nargin==4
    flag = 0;
end

N = size(G,1);
vertices = [];
k = 1;
for i = 1:N
    for j = i+1:N
        A = [G(i,:); G(j,:)];
        b = [g(i); g(j)];

        if rank(A) == 2
            vertex = A\b;
            if all(G*vertex<=(g+1e-3))
                vertices = [vertices;vertex'];
                k = k+1;
            end
        end
    end
end

% if size(vertices,1) ~= N
%     error('The htperplane plot error!')
% end

K = convhull(vertices);
vertices = vertices(K,:);

% 绘制矩形区域
if flag == 1
    figure;
    % grid on;
    xlabel('x轴');
    ylabel('y轴');
    title('四个超平面组成的矩形区域');
else
    hold on
end
fill(vertices(:, 1), vertices(:, 2), color,'FaceAlpha', FaceAlpha);
axis equal
grid on;
% axis([-2 2 -2 2]); % 设置坐标轴范围，适当调整以便更清楚显示正方形

end