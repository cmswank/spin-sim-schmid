Tflush=4;
Sx=zeros(size(z,2),1);
Sy=zeros(size(z,2),1);
Sz=zeros(size(z,2),1);
T=zeros(size(z,2),1);
%exind=find(t>Tflush,1,'first');
for i = 1:size(z,2)
exind=find(z(:,i)>1.9,1,'first');
    if isempty(exind)
        T(i)=Tflush+1;
        continue
    end
Sx(i)=sx(exind,i);
Sy(i)=sy(exind,i);
Sz(i)=sz(exind,i);
% Sx(i)=sx(end,i);
% Sy(i)=sy(end,i);
% Sz(i)=sz(end,i);

T(i)=t(exind);
end

%What is the polarization if we close the valve after Tflush seconds? 
disp(['Flush time: ',num2str(Tflush),' s, Polarization: ',num2str(mean(Sx(T<Tflush))),', number fraction: ',num2str(sum(T<Tflush)/size(sx,2))]);
%disp(['Flush time: ',num2str(Tflush),' s, Polarization: ',num2str(mean(Sx))]);

