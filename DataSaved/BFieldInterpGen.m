xpoints=11;
ypoints=11;
zpoints=11;
Bdata=zeros(xpoints*ypoints*zpoints,1);

ind=1;
 for i = 1:xpoints
    for ii=1:ypoints
        for iii=1:zpoints
            Bdata(ind)=3E-6; %Homogenous field for test. 
            ind=ind+1;
        end
    end
 end
%GOOD LUCK!!!


fileID = fopen('/data1/cmswank/spin-sim-xliu/BField/Bdatax.dat','w');
fwrite(fileID,Bdata,'double');
fclose(fileID);




