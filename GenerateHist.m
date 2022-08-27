function hist = GenerateHist(Bias_rest, LogG_rest, GateV_start, GateV_end, logG_start, logG_end)
%生成用于画图的数据

n_bins = 300;
% GateV_start = -0.6;
% GateV_end = 1;
% 
% logG_start = -8;
% logG_end = -4;

% HistData = histogram2(Bias_rest{1}, LogG_rest{1},'DisplayStyle','tile','ShowEmptyBins','on');
% histogram2(Bias_rest{1}, LogG_rest{1},'DisplayStyle','tile','ShowEmptyBins','on');
% HistData = histcounts2

Xedges = linspace(GateV_start,GateV_end,n_bins+1);
Yedges = linspace(logG_start,logG_end,n_bins+1);   %creat same length with PSD and Gavg

trace = length(Bias_rest);
hist = zeros(n_bins);
if length(Bias_rest) ~= 0
    for i=1:trace
        histcount0 = histcounts2(Bias_rest{i}, LogG_rest{i}, Xedges, Yedges);
        hist_temp = histcount0';
        hist = hist + hist_temp;
    end
end
