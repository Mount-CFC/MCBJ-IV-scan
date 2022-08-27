function [logG_open] = cutconductancestep(data_all,additionallength1)


%%%计算%%%%%%
%%%Parameter setting  @@@@@@@（STM-LGP）@@@@@@@@@@@
V0=0.1;
G0=12.9;          %（为1/G0）
zero_set=-0.3;
rate=0.0001;
add1=additionallength1 / 2;
% %@@@@@@@@@@@@@@@@@@@@ 拟合参数
offset=0;
% a1= -9.1137;  b1= -27.646;  c1= -1.1614E-11;   d1= 4.1597E-12;
% a2=9.2183;   b2= -27.8018;  c2= 1.1899E-11;  d2= -1.4714E-13;

a1 = -1.6383; b1 = -14.809 ;a2=1.638; b2 = -14.78; 
%%%VtoI
data_all=data_all-offset;   %correcting voltage by offset
% for i=1:size(data_all,2)
%     if data_all(i)<0
%         Cur0(i)=c1*data_all(i)+exp(data_all(i)*a1+b1)+d1;
%     else
%         Cur0(i)=c2*data_all(i)+exp(data_all(i)*a2+b2)+d2;
%     end
% end
% log(I) = a1 * v + b
point = size(data_all,2)     %仅查询第二个维度的长度
% 将放大器的电压转换为电流
for i=1:point
    if data_all(i)<0
        Cur0(i)=10^(data_all(i)*a1+b1);
    else
        Cur0(i)=10^(data_all(i)*a2+b2);
    end
end

% for i=1:size(Cur0,2)        % Transfer I to G
%     logG(i)=log10(abs(G0/((V0/1000)/Cur0(i)-0.1)));
% end
% Transfer I to G
for i=1:size(Cur0,2)        
    logG(i)=log10(abs(G0/((V0/1000)/Cur0(i))));
end

%%%截取%%%%
for i=1:size(data_all,2)-4  % Take four sampling points as
    mean_cur(i)=mean(Cur0(i:i+4));
end
k=1;
for i=1:size(data_all,2)-5
    if(mean_cur(i)<5e-5&&mean_cur(i+1)>5e-5) %find the close connecting point
        Close(k)=i;
        k=k+1;
    end
end

k=1;
for i=1:size(data_all,2)-5
    if(mean_cur(i)>5e-5&&mean_cur(i+1)<5e-5) %find the open snaping point
        Open(k)=i;
        k=k+1;
    end
end
for i=1:size(data_all,2)-20
    std_cur(i)=std(Cur0(i:i+20)); %calculate standard derivation between 20 points
end

if length(Open)~=length(Close)   % make the length of Open and Close be equal
    if length(Open)>length(Close)
        Open((length(Close)+1):end)=[];
    else
        Close((length(Open)+1):end)=[];
    end
end

if(Open(1)<Close(1)) %judge the type of the inital trace
    for i=1:length(Open) %for the open initial

        Wait(i)=max(find(std_cur(Open(i):Close(i))==min(std_cur(Open(i):Close(i))))); %find the background strating point
        logG_openA{i}=logG(Open(i):(Wait(i)+Open(i)));%cut the open trace and save in the struct logG.openA
        length_open(i)=floor(Wait(i)/1000)*1000;     %why? it's just the length of Wait(i)
    end
else
    
    for i=1:length(Open)-1

        Wait(i)=max(find(std_cur(Open(i):Close(i+1))==min(std_cur(Open(i):Close(i+1)))));
        logG_openA{i}=logG(Open(i):(Wait(i)+Open(i)));
        length_open(i)=floor(Wait(i)/1000)*1000;
    end
end

%%%筛选有addtionallength
z=1;
k=1;
if(Open(1)<Close(1))
    for i=1:length(Open)

        if(length_open(i)+add1>additionallength1)           % extending the length of open trace
            logG_open{k}=logG(Open(i):(Wait(i)+Open(i))+(additionallength1-length_open(i)));
            k=k+1;
        end
        
        
    end
else
    for i=1:length(Open)-1

        if(length_open(i)+add1>additionallength1)
            logG_open{k}=logG(Open(i):(Wait(i)+Open(i))+(additionallength1-length_open(i)));
            k=k+1;
        end
    end
end






