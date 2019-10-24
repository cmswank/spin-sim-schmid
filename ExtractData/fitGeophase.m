function [cfun,gof] = fitGeophase(phi,t)
%linear fitter,, nothing special. S must be normalized. 

fo = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',-0.5,...
               'Upper',0.5,...
               'StartPoint',1e-4);
ft = fittype('a*x','options',fo);

[cfun,gof]=fit(t,phi,ft);

end

