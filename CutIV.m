function [ForwardTraceBias,ForwardTraceCurrent,ForwardTraceLogG,ReverseTraceBias,ReverseTraceCurrent,ReverseTraceLogG] = CutIV(bias, current, logG, voltage)
%voltge: scan voltage range, usually -1 to +1

points = length(bias);

%截取开始的点,delete points before check points
for i = 1:points
    if bias(i) == voltage || bias(i) == -voltage
        StartPoint = i;
        break
    else
        StartPoint = 0;
    end
end

%截取结束的点,delete points after check points
for j = points:-1:1
    if bias(j) == voltage || bias(j) == -voltage
        EndPoint = j;
        break
    else
        EndPoint = 0;
    end
end

%截取每组完整扫IV的片段,complete IV curves
biasNew = bias(StartPoint : EndPoint);
currentNew = current(StartPoint : EndPoint);
logGNew = logG(StartPoint : EndPoint);

%找出所有正扫、反扫的点的索引
DiffBias = diff(biasNew);
ForwardIndex = find(DiffBias > 0);
ReverseIndex = find(DiffBias < 0);
% ForwardBias = biasNew(ForwardIndex);
% ForwardCurrent = currentNew(ForwardIndex);
% ReverseBias = biasNew(ReverseIndex);
% ReverseCurrent = currentNew(ReverseIndex);

%找出所有正扫、反扫的曲线的开始曲线的索引
DiffIndex_f = diff(ForwardIndex);
DiffIndex_r = diff(ReverseIndex);

ForwardTraceIndex = [0 find(DiffIndex_f>1000) length(DiffIndex_f)];
ReverseTraceIndex = [0 find(DiffIndex_r>1000) length(DiffIndex_r)];

%生成每条曲线构成的元胞数组
%正扫
for m = 1:length(ForwardTraceIndex)-1
    
    ForwardTraceBias{m} = biasNew(ForwardIndex(ForwardTraceIndex(m)+1) : ForwardIndex(ForwardTraceIndex(m+1)));
    ForwardTraceCurrent{m} = log10(abs(currentNew(ForwardIndex(ForwardTraceIndex(m)+1) : ForwardIndex(ForwardTraceIndex(m+1)))) .* 1e6);
    ForwardTraceLogG{m} = logGNew(ForwardIndex(ForwardTraceIndex(m)+1) : ForwardIndex(ForwardTraceIndex(m+1)));
%     if m > 30 && m < 40
%         figure(m+10)
% %         plot(ForwardTraceBias{m}, ForwardTraceCurrent{m})
%         plot(ForwardTraceCurrent{m})
%         title(num2str(m))
%     end
end
%反扫
for n = 1:length(ReverseTraceIndex)-1
%     figure(n)
    ReverseTraceBias{n} = biasNew(ReverseIndex(ReverseTraceIndex(n)+1) : ReverseIndex(ReverseTraceIndex(n+1)));
    ReverseTraceCurrent{n} = log10(abs(currentNew(ReverseIndex(ReverseTraceIndex(n)+1) : ReverseIndex(ReverseTraceIndex(n+1)))) .* 1e6);
    ReverseTraceLogG{n} = logGNew(ReverseIndex(ReverseTraceIndex(n)+1) : ReverseIndex(ReverseTraceIndex(n+1)));
%     plot(ReverseTraceCurrent{n})
end




