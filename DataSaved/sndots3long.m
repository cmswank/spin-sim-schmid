function [spindot] = sndots3long(sx1,sy1,sz1,sx2,sy2,sz2)
%GETTHETA Summary of this function goes here
%   Detailed explanation goes here
    spindot=sx1.*sx2+sy1.*sy2+sz1.*sz2;
    
    %theta=acos(spindot);
end

