function w_rf_add = ModulationFunc( rad, n_pow, scale )

%w_rf_add = (cos(rad)+scale*cos(rad*2)).^n_pow;
x = mod(rad,2*pi);
w_rf_add = 1/scale*exp(-n_pow*(x-pi/2).^2)-scale*exp(-n_pow*(x-3*pi/2).^2);

end

