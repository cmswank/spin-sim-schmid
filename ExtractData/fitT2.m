function [cfun,gof] = fitT2(t,S)
%linear fitter,, nothing special. S must be normalized. 
S=abs(S);

fo = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',0.001,...
               'Upper',1E9,...
               'StartPoint',1000);
ft = fittype('exp(-x/a)','options',fo);


[cfun,gof]=fit(t,S,ft);

end

