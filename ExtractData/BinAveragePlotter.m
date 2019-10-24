bindata00=zeros(1000,1);
for i=1:1000
    bindata00(i)=sum(binAverage(Signal4900(:,7)-Signal4900(:,9),i));
end