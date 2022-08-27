function [biasSelected, currentSelected]=SelectIV(bias, current)
% 尝试筛选出有起跳的IV
%paras：
% bias：电压曲线构成的元胞数组
% current：lg (nA)电流曲线构成的元胞数组

biasSelected = {};
currentSelected = {};
k=1;
for i=1:length(bias)
    if mean(current{i}(end-400:end)) >= 1
        biasSelected{k} = bias{i};
        currentSelected{k} = current{i};
        k = k+1;
    end
end