function [b] = BIN(dataset1,dataset2,nb)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
X=linspace(min(dataset1(dataset1>0)),max(dataset1(dataset1>0)),nb+1);

for i=1:nb
    in=find(and(dataset1>=X(i),dataset1<X(i+1)));
    x(i)=mean([X(i),X(i+1)]);
    mn(i)=mean(dataset2(in));
    med(i)=median(dataset2(in)');
    sd(i)=std(dataset2(in));
    IQ(i)=iqr(dataset2(in));
    % mx(i)=max(dataset2(in));
    % MIN(i)=min(dataset2(in));


end

b(:,1)=x';
b(:,2)=mn';
b(:,3)=med';
b(:,4)=sd;
% b(:,5)=MIN;
% b(:,6)=mx';
b(:,7)=IQ';


end