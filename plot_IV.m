function [hist] = plot_IV(x,y,xlow,xhigh,ylow,yhigh,xbins,ybins)
%生成x，y的二维统计直方图的通用函数
%x:x轴组成的元胞数组，每一个元胞为一条x轴的曲线
%y:y轴组成的元胞数组，每一个元胞为一条y轴的曲线
%return----画出一副hist2D

% 这里精修图片的坐标轴大小
Xstart = -1;
Xend = 1;
Ystart = ylow;
Yend = yhigh;


Xedges = linspace(xlow,xhigh,xbins+1);
Yedges = linspace(ylow,yhigh,ybins+1);   %creat same length with PSD and Gavg

trace = length(x);
hist = zeros(xbins,ybins);
if length(x) ~= 0
    for i=1:trace
        histcount0 = histcounts2(x{i}, y{i}, Xedges, Yedges);
        hist_temp = histcount0';
        hist = hist + hist_temp;
    end
end

Xaxis = linspace(xlow,xhigh,xbins);
Yaxis = linspace(ylow,yhigh,ybins);
imagesc(Xaxis,Yaxis, hist);
% imagesc(Xedges,Yedges, histcount);
set(gca,'YDir','normal')
set(gca,'tickdir','out')
load MyColormapRandB.mat;
% load Mycolormaps.mat;
colormap(mycmap);
colorbar;
% ylabel('Conductance / log (\itG/\itG\rm_0)', 'Interpreter', 'tex','FontSize',15)
xlabel('Voltage / V', 'Interpreter', 'tex','FontSize',15)
ylabel({'Current / lg (nA)'},'Interpreter','tex','FontSize',15)
xlim([Xstart, Xend]);
ylim([Ystart, Yend]);