clc
tic
clear all
close all
load traces_open.mat;
%rawtraceG=cell(1,length(traces_open));


highend=-2.5;     %@@@@@@@@@@@@ High End @@@@@@@@@@@@
lowend=-4;        %@@@@@@@@@@@@ Low End @@@@@@@@@@@@@@@


traces_open = logG_open;
for i=1:length(traces_open) % isolate G and L to rawtraceL and rawtraceG
    rawtraceG{i}=traces_open{i};
end

for i=1:length(traces_open)
    for k=1:length(traces_open{i})
        if rawtraceG{i}(k) < highend          %@@@@@@@@@@@@ High End @@@@@@@@@@@@
            cut_start(i) = k;
            break
        end
    end
end

traces_open(:)=[]; % release RAM

for i=1:length(rawtraceG)
    for k=1:length(rawtraceG{i})
        if rawtraceG{i}(k)<lowend           %@@@@@@@@@@@@ Low End @@@@@@@@@@@@@@@
            cut_end(i)=k;
            break
        end
    end
end

Lcycle=min(length(cut_end),length(cut_start));

for i=1:Lcycle % calculate G_AVG and extracting analyzing data

    if cut_start(i)>=cut_end(i);
        raw_for_flk{i}=[];
    else
        raw_for_flk{i}=rawtraceG{i}(cut_start(i):cut_end(i));
    end
    
    fprintf('Percentage: %03u\n',floor(i/length(rawtraceG)*100));
end

save raw_for_flk.mat raw_for_flk;

figure(2)
for i=1:10
    
    subplot(2,5,i);
    n=unidrnd(length(raw_for_flk))
    plot(raw_for_flk{n});
end

toc


    