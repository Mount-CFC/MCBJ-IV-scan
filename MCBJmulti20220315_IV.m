%适用于IV
% 一张图生成多张子图的logG-t图用于粗筛,记录AO1 bias，AI0 sampling Voltage, current, Conductance LogG
%同时画logG与current曲线
%画log I 曲线

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

for n = 1:num_files
    struc=TDMS_readTDMSFile(filename1{n});
    data_bias=struc.data{1,3};              %第一行第3列，提取Bias
    data_samplingV = struc.data{1,4};       %第一行第4列，提取sampling voltage
    data_Cur = struc.data{1,5};             % 第一行第5列，提取current
    data_logG = struc.data{1,6};            % 第一行第6列，提取log (G/G0)
%     转换为电流
%     Cur = zeros(size(data_samplingV));
%     Cur(data_samplingV >= 0) = 10.^(data_samplingV(data_samplingV >= 0)*a2 + b2);
%     Cur(data_samplingV < 0) = 10.^(data_samplingV(data_samplingV < 0)*a1 + b1);
%     转换为电导
%     logG = log10(Cur ./ data_bias / 77.6e-6);

    %Conductance plot
    figure(1)
    subplot(num_files,1,n)
%     https://blog.csdn.net/qingfengxd1/article/details/120122017  绘制双y轴
    
    x = 1:length(data_bias);
    yyaxis left
    %conductance
    plot(x, data_logG)
    title(filename1{n},'FontSize',5)
    ylabel('Conductance / log (\itG/\itG\rm_0)', 'Interpreter', 'tex','FontSize',5)
    xlabel({'Sampling points / 100us per point'},'Interpreter','tex','FontSize',5)
    ylim([-5.5 -2.5])
    
    yyaxis right
    %bias
    plot2 = plot(x, data_bias,'LineWidth', 2, 'LineStyle','-');
    plot2.Color(4) = 0.6;       %调节Color(4)这个参数可以设置不同的透明度
    hold on
    line([0 length(data_bias)], [0 0], 'linestyle', '--', 'Color', 'r', 'LineWidth', 2)
    ylabel('Bias / V', 'FontSize' ,5)
    ylim([-2 2])
    
    %Current Plot
    figure(2)  
    subplot(num_files,1,n)
%     https://blog.csdn.net/qingfengxd1/article/details/120122017  绘制双y轴
    
    x = 1:length(data_bias);
    yyaxis left
%     current
    plot(x, data_Cur)
    title(filename1{n},'FontSize',5)
    ylabel('Current / mA','FontSize',5)
%     ylabel('Current / lg(nA)','FontSize',5)

    xlabel({'Sampling points / 100us per point'},'Interpreter','tex','FontSize',5)             
%     ylim([-2e-5 2e-5])
    
    yyaxis right
    %bias
    plot3 = plot(x, data_bias,'LineWidth', 2, 'LineStyle','-');
    plot3.Color(4) = 0.6;       %调节Color(4)这个参数可以设置不同的透明度
    hold on
    line([0 length(data_bias)], [0 0], 'linestyle', '--', 'Color', 'r', 'LineWidth', 2)
    ylabel('Bias / V', 'FontSize' ,5)
    ylim([-10 10])
    
    %Lg (Current) Plot
    figure(3)
    subplot(num_files,1,n)
%     https://blog.csdn.net/qingfengxd1/article/details/120122017  绘制双y轴
    
    x = 1:length(data_bias);
    yyaxis left
%     current
%     plot(x, data_Cur)
    %以nA为单位，取对数得到电流值
    log_cur = log10(abs(data_Cur) .* 1e6);
    plot(x, log_cur);
    title(filename1{n},'FontSize',5)
%     ylabel('Current / mA','FontSize',5)
    ylabel('I / lg(nA)','FontSize',5)

    xlabel({'Sampling points / 100us per point'},'Interpreter','tex','FontSize',5)             
%     ylim([-2e-5 2e-5])
    
    yyaxis right
    %bias
    plot3 = plot(x, data_bias,'LineWidth', 2, 'LineStyle','-');
    plot3.Color(4) = 0.6;       %调节Color(4)这个参数可以设置不同的透明度
    hold on
    line([0 length(data_bias)], [0 0], 'linestyle', '--', 'Color', 'r', 'LineWidth', 2)
    ylabel('Bias / V', 'FontSize' ,5)
    ylim([-10 10])
    
    %     clear test data_s
end




toc