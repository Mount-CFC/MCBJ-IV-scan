clc
clear
close all
tic


fs=10000; % @@@ Sampling Frequency @@@10kHz采样�?
ft=1/fs;
traceL = 10000  %%@@@@@@@@@@ traces length (how many points per trace)@@@@@@@@@@@@@
lowG=-2.8;
highG=-2.5;
lowPSD=-6;
highPSD=-5;
n_bins=50;         %热力图的分辨率，做histcount2的格子数�?
PlotSelect = 0     %If it is True, polt the raw .TDMS files

[filename,filepath]=uigetfile('*.tdms','Select data files','MultiSelect','on');
if iscell(filename)
    filename1=filename;
else
    filename1{1}=filename;
end
num_files = length(filename1)

%%% plot long traces
if PlotSelect == 1
    figure(20)
    if num_files == 1
        test=TDMS_readTDMSFile(filename1{1});
        data_s=test.data{1,4};
        plot(data_s)
        title(filename1{1},'FontSize',5)
        ylabel('Conductance / log (\itG/\itG\rm_0)', 'Interpreter', 'tex','FontSize',5)
        xlabel({'Sampling points / 50us per point'},'Interpreter','tex','FontSize',5)
    else  
        for n = 1:num_files
            test=TDMS_readTDMSFile(filename1{n});
            data_s=test.data{1,4}; 
            subplot(num_files,1,n)
            plot(data_s)
            title(filename1{n},'FontSize',5)
            ylabel('Conductance / log (\itG/\itG\rm_0)', 'Interpreter', 'tex','FontSize',5)
            xlabel({'Sampling points / 50us per point'},'Interpreter','tex','FontSize',5)
        %     clear test data_s

        end
    end
    saveas(gcf,'1_TraceBeforeCut.fig')
end

%% 获得原始数据
data_s = [];
if num_files == 1
    test=TDMS_readTDMSFile(filename1{1});
    data_s=test.data{1,4};
else  
    for n = 1:num_files
        test=TDMS_readTDMSFile(filename1{n});
        data_temp=test.data{1,4}; %第一行第4列，提取Conductance
        data_s = [data_s,data_temp];
        disp(['loading ' filename1{n} '...']); % Present the file name
    %     clear test data_s

    end
end
%%

for i=1:(length(data_s)/traceL)
    raw_for_flk{i} = data_s(:, (traceL*(i-1)+1):traceL*i);
end
%% 生成单条原始数据的预�?
figure(1)
for i=1:10
    
    subplot(5,2,i);
    n=unidrnd(length(raw_for_flk));
    plot(raw_for_flk{n});
    title(n)
%     ylim([lowG highG])
    ylim([-3.6 -2.4])
end
saveas(gcf,'2_RawDataForPSD.fig')

%% @@@@@@@@@@@@@@@@@@@@@@@@@@@@@ FFT and calculate PSD
disp(['Data processing...'])
for k=1:length(raw_for_flk)
    data_L=2^nextpow2(length(raw_for_flk{k})); %   Determine Number of points (next power of 2)
    yt=transpose(raw_for_flk{k});  % Transpose the column vector to row vector
    yt=power(10,yt.*1);
    Fyt=fft(yt,data_L)/data_L*2; % FFT amplitude
    f=fs/data_L*(0:1:data_L-1); % Frequency range
    F=fs/data_L;
    A_Fyt=abs(Fyt);
    SQ_Fyt=A_Fyt.^2; % square
    
    for z=0:data_L              %  intergration strating range
        if fs/data_L*z > 100
            sp=z;
            break
        end
    end
    
    for z=sp:data_L              %  intergration end range
        if fs/data_L*z > 1000
            ep=z;
            break
        end
    end
    %int_range=(100:fs/data_L:1000);
    PSD0(k)=trapz(f(sp:ep),SQ_Fyt(sp:ep));            % @@@ Trapezoidal integration
    % @@@ correction integration  100 and 1000 are the intergation range
    if f(sp)~=100
        PSD0(k)=PSD0(k)+0.5*(((SQ_Fyt(ep+1)-SQ_Fyt(ep))/F*(1000-f(ep))+SQ_Fyt(ep))+SQ_Fyt(ep))*(1000-f(ep))-0.5*(((SQ_Fyt(sp+1)-SQ_Fyt(sp))/F*(100-f(sp))+SQ_Fyt(sp))+SQ_Fyt(sp))*(100-f(sp));
    end
    MeanG(k)=log10(mean(yt));
    PSD(k)=log10(PSD0(k)/(10^MeanG(k)));
    % PSD(k)=PSD(k)/(10^MeanG(k));
    %figure(k);
    % plot(f(1:data_L/2),SQ_Fyt(1:data_L/2));
    % xlim([90,1100]);
end

%maxpsd= max(PSD);
%minpsd= min(PSD);
n=1;
meanlogG=mean(MeanG);
%% 绘制PSD热力�?

% [histcount,Xedges,Yedges] = histcounts2(MeanG,PSD,[n_bins,n_bins]);


% x = randn(1000,1);
% y = randn(1000,1);
% Xedges = -5:5;
% Yedges = [-5 -4 -2 -1 -0.5 0 0.5 1 2 4 5];
% N = histcounts2(x,y,Xedges,Yedges)


figure(2)
Xedges = linspace(lowG,highG,n_bins+1);
Yedges = linspace(lowPSD,highPSD,n_bins+1);   %格子数目比hist数据多一�?
% h = histogram2(MeanG,PSD,Xedges,Yedges,'DisplayStyle','tile','ShowEmptyBins','on','LineStyle','none')
% histcounts = h.Values'

% histcount = histcounts2(MeanG,PSD,Xedges,Yedges)
% Xaxis = linspace(lowG,highG,n_bins);
% Yaxis = linspace(lowPSD,highPSD,n_bins);


% imagesc(Xaxis,Yaxis, histcount);
% imagesc(Xedges,Yedges, histcount);
% set(gca,'YDir','normal')
% set(gca,'tickdir','out')
% load MyColormapRandB.mat;
% colormap(mycmap);
% colorbar;
% xlabel('Conductance / log (\itG/\itG\rm_0)', 'Interpreter', 'tex') 
% ylabel('Noise Power / log (\itG/\itG\rm_0)', 'Interpreter','tex')

%%%%% test
% figure(100)
% Xaxis = linspace(lowG,highG,n_bins);
% Yaxis = linspace(lowPSD,highPSD,n_bins);
% 
% 
% imagesc(Xaxis,Yaxis, histcounts);
% % imagesc(Xedges,Yedges, histcount);
% set(gca,'YDir','normal')
% set(gca,'tickdir','out')
% load MyColormapRandB.mat;
% colormap(mycmap);
% colorbar;
% xlabel('Conductance / log (\itG/\itG\rm_0)', 'Interpreter', 'tex') 
% ylabel('Noise Power / log (\itG/\itG\rm_0)', 'Interpreter','tex')


% figure(101)
histcount0 = histcounts2(MeanG,PSD,Xedges,Yedges)
Xaxis = linspace(lowG,highG,n_bins);
Yaxis = linspace(lowPSD,highPSD,n_bins);
histcount = histcount0';

imagesc(Xaxis,Yaxis, histcount);
% imagesc(Xedges,Yedges, histcount);
set(gca,'YDir','normal')
set(gca,'tickdir','out')
load MyColormapRandB.mat;
colormap(mycmap);
colorbar;
xlabel('Conductance / log (\itG/\itG\rm_0)', 'Interpreter', 'tex','FontSize',20,'FontName','Arial') 
ylabel('Noise Power / \itG \rm(log \itG\rm_0)', 'Interpreter','tex','FontSize',20,'FontName','Arial')
set(gca,'FontSize',15,'LineWidth',1.5,'FontName','Arial')

%% 椭圆拟合
hold on;
figure(2)
[fitresult0, gof] = createFit(Xaxis, Yaxis, histcount);

xRangeStart=lowG;
xRangeEnd=highG;
yRangeStart=lowPSD;
yRangeEnd=highPSD;
    
    
a =      fitresult0.a ;
b =      fitresult0.b;
c =      fitresult0.c  ;
d =      fitresult0.d  ;
e =      fitresult0.e  ;
f =      fitresult0.f  ;

meanshift=(1-n)*meanlogG;
[x,y]=meshgrid(xRangeStart:0.01:xRangeEnd,(yRangeStart+meanshift):0.01:(yRangeEnd+meanshift));
z = exp(a*x.^2 + x.*y.*b + x.*c + y.*d + e*y.^2 + f);
contour(x,y,z,4,'-k','linewidth',1.5);
saveas(gcf,'3_PSD.fig')
    
    
    
    
%%
% 可能的改进：
% 筛�?�出好的电导区间(通过抖动、不抖动求出点数)
% 确定好椭圆拟合的的range，试试用PSD统一，设置好psd high low


%% 确定n的�??
% clear histcount histcount0
% lowG=-2.8;
% highG=-2.4;
lowPSDForN=-5.5;
highPSDForN=-4.7;
% Xedges = linspace(lowG,highG,n_bins+1);
Yedges = linspace(lowPSDForN,highPSDForN,n_bins+1);   %格子数目比hist数据多一�?
% Xaxis = linspace(lowG,highG,n_bins);
Yaxis = linspace(lowPSDForN,highPSDForN,n_bins);

for n=1.1:0.1:2
    meanshift=(1-n)*meanlogG;
    for k=1:length(PSD)
        PSD1(k)=log10(PSD0(k)/(10^MeanG(k))^n);
    end

%     PSD_data1= [transpose(MeanG),transpose(PSD1)];
%     [PSD_2D11,X11,Y11]=plotPSD(PSD_data1,meanshift+highPSD,meanshift+lowPSD,highG,lowG,n_bins); %(data,PSD_max,PSD_min,G_max,G_min)
%     [histcount,Xedges,Yedges] = histcounts2(MeanG,PSD1,[n_bins,n_bins]);

    histcount0 = histcounts2(MeanG,PSD1,Xedges,Yedges);
    histcount = histcount0';
    [fitresult0, gof] = createFit(Xaxis, Yaxis, histcount);

    figure(round(2+10*(n-1)));

    imagesc(Xaxis, Yaxis, histcount);
    set(gca,'YDir','normal')
    set(gca,'tickdir','out')
    load MyColormapRandB.mat;
    colormap(mycmap);
    colorbar;
    xlabel('Conductance / log (\itG/\itG\rm_0)', 'Interpreter', 'tex','FontSize',20,'FontName','Arial') 
    ylabel(['Noise Power / \itG \rm(log \itG\rm_0^{n})  n=' num2str(n)], 'Interpreter','tex','FontSize',20,'FontName','Arial')
    set(gca,'FontSize',15,'LineWidth',1.5,'FontName','Arial')


hold on;
% xRangeStart=Xedges(1);
% xRangeEnd=Xedges(n_bins);
% yRangeStart=Yedges(1);
% yRangeEnd=Yedges(n_bins);
xRangeStart=lowG;
xRangeEnd=highG;
yRangeStart=lowPSD;
yRangeEnd=highPSD;


           a =      fitresult0.a ;
           b =      fitresult0.b;
           c =      fitresult0.c  ;
           d =      fitresult0.d  ;
           e =      fitresult0.e  ;
           f =      fitresult0.f  ;

        meanshift=(1-n)*meanlogG;
        [x,y]=meshgrid(xRangeStart:0.01:xRangeEnd,(yRangeStart+meanshift):0.01:(yRangeEnd+meanshift));
        z = exp(a*x.^2 + x.*y.*b + x.*c + y.*d + e*y.^2 + f);
        contour(x,y,z,5,'-k','linewidth',1.5);
end


%@@@@@@@@@@@@@@@@@@@  ����PSD��meanG�����ϵ�����ҵ���С�ĳ߶Ȼ�����


coe=[];
k=1;
for i = 0.5:0.1:2.5
   
   PSD_temp=log10(PSD0./(10.^MeanG).^i);
   coe_tmp=corrcoef(PSD_temp,MeanG);
   coe(k,1:2)=[i, coe_tmp(1,2)];
   

   k=k+1;
   
end
min_corrcoef=coe(abs(coe(:,2))==min(abs(coe(:,2))));

disp(['The scaling exponent zero correlation between the normalized ficker noise power and mean conductance: ',num2str(min_corrcoef)]); 





toc