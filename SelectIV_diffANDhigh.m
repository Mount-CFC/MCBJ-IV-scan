function [biasSelected2, currentSelected2]=SelectIV_diffANDhigh(bias, current,low, high)
% 使用diff+大于0.8的电压有一个高电流阈值，尝试筛选出有起跳的IV
%paras：
% bias：电压曲线构成的元胞数组
% current：lg (nA)电流曲线构成的元胞数组
%low:Average conductance should LOWER than this value between -0.3~-0.2V
% high:Average conductance should HIGHER than this value between -0.9~-0.8V


biasSelected = {};
currentSelected = {};
biasSelected2 = {};
currentSelected2 = {};

%正偏压下有起跳
% k=1;
% for i=1:length(bias)
% %     range = current{i}(bias{i} >= 0.1)
%     if max(diff(current{i}(bias{i} <= -0.1))) >= 0.01
%         biasSelected{k} = bias{i};
%         currentSelected{k} = current{i};
%         k = k+1;
%     end
% end


%低偏压下有低导,lowcond1\highcond1设置低导筛选区间
lowcond1 = -0.3;
highcond1 = -0.2;
k=1;
for i=1:length(bias)
    CurReal = (10 .^ current{i}) .* 1e-9;
    logG = log10((CurReal ./ abs(bias{i})) ./ (77.6e-6));
%     range = current{i}(bias{i} >= 0.1)
    if mean(logG(bias{i} <= highcond1 & bias{i} >= lowcond1)) <= low
        biasSelected{k} = bias{i};
        currentSelected{k} = current{i};
        k = k+1;
    end
    clear CurReal logG
end

% 大于0.9V有高电导，lowcond2\highcond2设置高导筛选区间
lowcond2 = -0.9;
highcond2 = -0.8;
l=1;
for i=1:length(biasSelected)
    CurReal = (10 .^ currentSelected{i}) .* 1e-9;
    logG = log10((CurReal ./ abs(biasSelected{i})) ./ (77.6e-6));
%     range = current{i}(bias{i} >= 0.1)
    if mean(logG(biasSelected{i} >= lowcond2 & biasSelected{i} <= highcond2)) >= high
        biasSelected2{l} = biasSelected{i};
        currentSelected2{l} = currentSelected{i};
        l = l+1;
    end
    clear CurReal logG
end