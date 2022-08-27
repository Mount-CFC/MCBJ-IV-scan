% 计算扫速,默认扫描范围为-1到+1
us = 100; %points interval

eq1 = find(abs(data_bias) == 1);
diff1 = diff(eq1);
FindPoint1 = diff1(diff1 ~= 1);
Points = mean(FindPoint1);
Rate = 1/(Points*us*1e-6/2);
fprintf('IV rate:%f V/s\n', Rate);