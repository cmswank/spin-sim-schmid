%%SimSimMatlabImportNew.m

Bin = 1001;

Event = 1000;

fileID = fopen('/home/cmswank/spin-sim-heatflush/ExtractData/OutputHeatflush001.dat');

A = fread(fileID, 'double');

B = reshape(A, 7, Bin, Event);

sx = squeeze(B(1,:,2:Event));

sy = squeeze(B(2,:,2:Event));

sz = squeeze(B(3,:,2:Event));

x =  squeeze(B(4,:,2:Event));

y =  squeeze(B(5,:,2:Event));

z = squeeze(B(6,:,2:Event));


tlarge = squeeze(B(7,:,2:Event)); 

t = squeeze(tlarge(:,1));

