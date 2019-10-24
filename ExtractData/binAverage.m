function [data] = binAverage(inputdata,binsize)
    data=zeros(ceil(length(inputdata)/binsize),1);
    for i = 1:floor(length(inputdata)/binsize)
    data(i)=mean(inputdata(1+(i-1)*binsize:i*binsize));
    end
    if mod(length(inputdata),binsize)~=0
        data(i+1)=mean(inputdata(floor(length(inputdata)/binsize):end));
    end
end

