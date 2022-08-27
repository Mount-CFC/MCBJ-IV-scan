% select one sigle data to obtain the hist plot
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


test=TDMS_readTDMSFile(filename1{1}); %obtain a struct of TDMS file
data_s=test.data{1,5}; %Conductance is in the 5th column
data_bias = test.data{1,3}; %bias is in the 3rd column
% histogram(data_s, 1000)  %all histogram
lower_time = 30000; 
upper_time = 38500;
figure
subplot(121)
x = lower_time:upper_time;
yyaxis left
%Conductance
plot(x, data_s(:, lower_time:upper_time))
title(filename1{1})
ylabel('Conductance / log (\itG/\itG\rm_0)', 'Interpreter', 'tex','FontSize',12)
% title({['abc','L Range:',num2str(a),'(nm)','~~',num2str(b),'(cm)'];['B Range:',num2str(c),'(cm)','~~',num2str(d),'(cm)']})
% https://www.cnblogs.com/AI-Algorithms/p/3731232.html  
% Dont forget: num2str
xlabel({['Sampling points'];['From ' num2str(lower_time)  ' to '  num2str(upper_time)]},'Interpreter','tex','FontSize',12)
yyaxis right
%Bias
plot2 = plot(x, data_bias(:, lower_time:upper_time), 'LineWidth', 2, 'LineStyle','-');
plot2.Color(4) = 0.6; %调节Color(4)这个参数可以设置不同的透明度
ylabel('Bias / V', 'FontSize' ,15)
ylim([0 0.6])



subplot(122)
histogram(data_s(:,lower_time:upper_time), 1000,'Orientation','horizontal')
xlabel('Conductance / log (\itG/\itG\rm_0)', 'Interpreter', 'tex','FontSize',12)
ylabel({'Counts'},'Interpreter','tex','FontSize',12)
toc