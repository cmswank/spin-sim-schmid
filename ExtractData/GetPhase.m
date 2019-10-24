function [ phase ] = GetPhase( time )

global fm;
global wrf;
global deltat1;
global deltat2;
global dw1;
global dw2;
global phi;

halves = floor(time*2.0*fm); % no. of complete half periods that have passed
halfT = 0.5/fm;           % half the period
phase = phi + wrf*halves*halfT; % phase of complete half periods + phi
if mod(halves,2) == 0
    phase = phase + halves/2 * (dw1*deltat1 + dw2*deltat2); % no. of complete periods * 
else
    phase = phase + ( (halves-1)/2 *(dw1*deltat1 + dw2*deltat2)+ dw1* deltat1 ); 
end
mod_t = time - halfT * halves;
phase = phase + mod_t * wrf;   % adding normal unmodulated phase 
if mod(halves,2) == 0
    if mod_t < deltat1
        phase = phase + dw1* mod_t;
    else
        phase = phase + dw1* deltat1;
    end
else
    if mod_t < deltat2
        phase = phase + dw2* mod_t;
    else
        phase = phase + dw2* deltat1;
    end
end

end
