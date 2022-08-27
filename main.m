clc
close all
clear all
tic

additionallength1=20000;%@@@@@@@@@@@@@@@ 调整该值与labview处理系统一致
%这个值与logG_open的每个cell的length有关联

[filename,filepath]=uigetfile('*.tdms','Select data files','MultiSelect','on');
if iscell(filename)
    filename1=filename;
else
    filename1{1}=filename;
end

num_file = length(filename1);
fileclass = '.tdms';
logG_open=[];
data_s=cell(1,100);
for i = 1 : num_file
    test=TDMS_readTDMSFile(filename1{i});
    data_s=test.data{1,3};
    logG_open_trans=cutconductancestep(data_s,additionallength1);
    logG_open=[logG_open,logG_open_trans]
%     clear test data_s logG_open_trans;
    save traces_open.mat logG_open;
    fprintf('Number: %03u\n',floor(i)); % Present the number    
end

figure(1)
for i=1:10
    
    subplot(2,5,i);
    n=unidrnd(length(logG_open));
    plot(logG_open{n});
end


toc
