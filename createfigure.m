function createfigure(X1, Y1, Y2)
%CREATEFIGURE1(X1, Y1, Y2)
%  X1:  x 数据的矢量
%  Y1:  y 数据的矢量
%  Y2:  y 数据的矢量

%  由 MATLAB 于 15-Oct-2016 13:14:45 自动生成

% 创建 figure
figure1 = figure;

% 创建 subplot
subplot1 = subplot(2,1,1,'Parent',figure1);
%% 取消以下行的注释以保留坐标轴的 X 范围
% xlim(subplot1,[0.001 0.004]);
%% 取消以下行的注释以保留坐标轴的 Y 范围
% ylim(subplot1,[38 38.6]);
box(subplot1,'on');
hold(subplot1,'on');

% 创建 ylabel
ylabel('PSNR (dB)');

% 创建 xlabel
xlabel('\lambda');

        xlim(subplot1,[1.5e-3 2.5e-3]);
        ylim(subplot1,[38.4 38.6]);
% 创建 plot
plot(X1,Y1,'Parent',subplot1,'MarkerFaceColor',[0 0 0],'MarkerSize',3,...
    'Marker','o',...
    'LineWidth',1,...
    'LineStyle','-.',...
    'Color',[0 0 1]);

% 创建 subplot
subplot2 = subplot(2,1,2,'Parent',figure1);
%% 取消以下行的注释以保留坐标轴的 X 范围
% xlim(subplot2,[0.001 0.004]);
%% 取消以下行的注释以保留坐标轴的 Y 范围
% ylim(subplot2,[0.96 0.97]);
box(subplot2,'on');
hold(subplot2,'on');

% 创建 ylabel
ylabel('SSIM');

% 创建 xlabel
xlabel('\lambda');

   xlim(subplot2,[1.5e-3 2.5e-3]);
        ylim(subplot2,[0.960 0.970]);
% 创建 plot
plot(X1,Y2,'Parent',subplot2,'MarkerFaceColor',[0 0 0],'MarkerSize',3,...
    'Marker','o',...
    'LineWidth',1,...
    'LineStyle','-.',...
    'Color',[1 0 0]);

