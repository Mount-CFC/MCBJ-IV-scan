function [biasSelected, currentSelected]=SelectIV_diff(bias, current)
% 使用diff函数尝试筛选出有起跳的IV
%paras：
% bias：电压曲线构成的元胞数组
% current：lg (nA)电流曲线构成的元胞数组

biasSelected = {};
currentSelected = {};
k=1;
for i=1:length(bias)
%     range = current{i}(bias{i} >= 0.1)
    if max(diff(current{i}(bias{i} >= 0.1))) >= 0.5
        biasSelected{k} = bias{i};
        currentSelected{k} = current{i};
        k = k+1;
    end
end