function [T2,T2upper,T2lower] = sparseT2fit(t,S73,Event,T2low,T2high,len,n)
%n is resolution
 T2found=false;
 T2upperfound=false;
 T2lowerfound=false;
% n=1E5;
% S73=S1;
% Event=1000;
% T2low=100;
% T2high=1E9;

if nargin < 6
    len=0.75;
end
if nargin<7
    n=1E5;
end


if nargin<3
    Event=1000;
end

if nargin<5
    T2high=1E7;
end


if nargin<4
    T2low=10;
end

for i=1:n+1
T2=T2high+(i-1)*(T2low-T2high)/n;
S=exp(-t/T2);
dS=S-S73;
  
  if any(dS(round(length(t)*len):end)<0)
    if T2found==true
      break
    else
        warning('lower T2 limit');
        T2=nan;
        break
    end
  else
      T2found=true;
      if i==n+1
       warning('increase T2 upper limit');
       T2=nan;
        break
      end  
      
  end
  
end


for i=1:n+1
T2upper=T2high+(i-1)*(T2low-T2high)/n;
S=exp(-t/T2upper);
S72=1-(1-1/sqrt(Event))*(1-S73);

dS=S-S72;
  
  if any(dS(round(length(t)*len):end)<0)
    if T2upperfound==true
      break
    else
        warning('lower T2 limit');
        T2upper=nan;
        break
    end
  else
      T2upperfound=true;
      if i==n+1 
       warning('increase T2 upper limit');
       T2upper=nan;
        break
      end  
      
  end
  
end

for i=1:n+1
T2lower=T2high+(i-1)*(T2low-T2high)/n;
S=exp(-t/T2lower);
S71=1-(1+1/sqrt(Event))*(1-S73);
dS=S-S71;
  
  if any(dS(round(length(t)*len):end)<0)
    if T2lowerfound==true
      break
    else
        warning('lower T2 limit');
        T2lower=nan;
        break
    end
  else
      T2lowerfound=true;
      if i==n+1
       warning('increase T2 upper limit');
       T2lower=nan;
        break
      end  
      
  end
  
end

end


