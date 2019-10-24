function [theta] = getmeanTheta(sx1,sy1,sz1,sx2,sy2,sz2)
%GETTHETA Summary of this function goes here
%   Detailed explanation goes here
    spindot=mean(sx1,2).*mean(sx2,2)+mean(sy1,2).*mean(sy2,2)...
        +mean(sz1,2).*mean(sz2,2);
    theta=acos(spindot);
end

