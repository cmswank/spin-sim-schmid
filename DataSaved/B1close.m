trials=1000;
wrf=6283.18530718;
 t=0:2E-5:1;
 tsh=0:2E-4:1;
 [A,B,C,D] = butter(10,[40/2500 120/2500],'stop');
 sos=ss2sos(A,B,C,D);
 
 for i = 1:trials 


 
  
end
%GOOD LUCK!!!







