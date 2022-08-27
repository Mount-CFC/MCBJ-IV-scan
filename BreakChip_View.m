%适用于芯片顶断过程的绘图
%记录motor position, Conductance LogG


clc
clear 
close all
tic



[filename,filepath]=uigetfile('*.tdms','Select data files','MultiSelect','on');
if iscell(filename)
    filename1=filename;
else 
    filename1{1}=filename;
end

num_files = length(filename1)
%%
%记录的AO1 bias，AI0 sampling voltage，current,  Conductance LogG
% a2 = 4.1106; b2 = -13.993 ;a1=-4.1242; b1 = -14.017; %MCBJ-Raman
data_logG = [];
data_bias = [];
for n = 1:num_files
    struc=TDMS_readTDMSFile(filename1{n});
    data_bias_temp= struc.data{1,3};              %第一行第3列，提取Bias
    data_logG_temp = struc.data{1,4};            % 第一行第6列，提取log (G/G0)
    data_bias = [data_bias data_bias_temp];
    data_logG = [data_logG data_logG_temp];
end


%Conductance plot
figure(1)
%     https://blog.csdn.net/qingfengxd1/article/details/120122017  绘制双y轴

x = 1:length(data_bias);
yyaxis left
%conductance
plot(x, data_logG)
title(filename1{n},'FontSize',5)
%     ylabel('Conductance / log (\itG/\itG\rm_0)', 'Interpreter', 'tex','FontSize',15)
ylabel('Log \itG', 'Interpreter', 'tex','FontSize',15)
xlabel({'Motor position'},'Interpreter','tex','FontSize',15)


yyaxis right
%bias
plot2 = plot(x, data_bias,'LineWidth', 2, 'LineStyle','-');
plot2.Color(4) = 0.6;       %调节Color(4)这个参数可以设置不同的透明度
ylabel('count', 'FontSize' ,15)

 
toc