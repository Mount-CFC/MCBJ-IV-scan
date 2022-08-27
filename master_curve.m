function [Xfitted,Yfitted] = master_curve(hist, xlow,xhigh,ylow,yhigh,xbins,ybins)
%拟合masercurve,xlow,xhigh,ylow,yhigh,xbins,ybins与plot_IV前面一致
% parms：hist，根据histgram2求出的histcount；
% Xaxis，根据linspace求出的从-scanrange到+scanrange除以bins的均匀点的行向量
% return：hist每一列的高斯分布峰中心的index

ft = fittype( 'gauss1' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';

Xaxis = linspace(xlow,xhigh,xbins);
Yaxis = linspace(ylow,yhigh,ybins);


for i=1:length(Xaxis)
    if all(hist(:,i) == 0)
        peakmiddle(i)=0;
    else
        try
            [fitresult] = fit((1:length(hist(:,i)))', hist(:,i), ft, opts ); 
            % Plot fit with data.
            peakmiddle(i)=fitresult.b1;
            % 高斯分布的中心点
        catch
%             如果不能拟合出guess，则取最大值；两个最大值取最大值的平均值
            selonecolumn = hist(:,i);
            columMAX = max(selonecolumn);
            peakmiddle(i) = mean(find(selonecolumn == columMAX));
        end
    end
end

z=1;
for i=1:length(peakmiddle)
    if peakmiddle(i)~=0
        Yfitted(z)=ylow + (yhigh-ylow)*peakmiddle(i)/length(Yaxis);
        Xfitted(z)=Xaxis(i);% 将等于0的扣掉，防止plot()的时候x轴y轴不一致
        z=z+1;
   
    end
end
% Xfitted = Xaxis;

% plot(Xfitted,Yfitted,'-b','linewidth',2);