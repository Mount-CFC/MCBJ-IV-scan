% Select one or more matrix file(s) by hist2d()
clc
tic
clear
close all
% load traces_open.mat;
%rawtraceG=cell(1,length(traces_open));
filename = uigetfile('*.mat','Select traces_open1,2,3... files','MultiSelect','on');
if iscell(filename)
    filename1=filename;
else
    filename1{1}=filename;
end
%%%%这些值与main程序保持一致,将值复制过来即可
n_bins = 300;
logG_start = -8;        %电导（y轴）的取值范围，低导
logG_end = -2;          %电导（y轴）的取值范围，高导
GateV_start = -1;       %gate电压（x轴）扫描的区间开始
GateV_end = 1;          %gate电压（x轴）扫描的区间结束

m_begin = -4.5;         %开始的10个点需要大于这个电导
m_end = -5.5;           %结束的最后10个点需要大于这个电导

% 这里精修图片的坐标轴大小
Xstart = -1;
Xend = 1;
Ystart = -8;
Yend = -2;


num_file = length(filename1);  %number of files selected
% logG_open = [];
histcount = zeros(n_bins);
for i = 1:num_file
    load(filename1{i})
    histcount = histcount + MatrixData;
    clear MatrixData
end
    



Xaxis = linspace(GateV_start,GateV_end,n_bins);
Yaxis = linspace(logG_start,logG_end,n_bins);
imagesc(Xaxis,Yaxis, histcount);
% imagesc(Xedges,Yedges, histcount);
set(gca,'YDir','normal')
set(gca,'tickdir','out')
load MyColormapRandB.mat;
colormap(mycmap);
colorbar;
ylabel('Conductance / log (\itG/\itG\rm_0)', 'Interpreter', 'tex','FontSize',15)
xlabel({'Potential / V \itvs. \rmAg/AgCl)'},'Interpreter','tex','FontSize',15)
xlim([Xstart, Xend]);
ylim([Ystart, Yend]);
% ylabel('Noise Power / \itG \rm(log \itG\rm_0)', 'Interpreter','tex','FontSize',20,'FontName','Arial')

