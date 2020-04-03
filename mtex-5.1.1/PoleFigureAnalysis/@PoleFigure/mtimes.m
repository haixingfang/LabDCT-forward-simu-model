function pf = mtimes(arg1,arg2)
% scaling of PoleFigures, implements pf1 * b and a * pf2
%
% overload the * operator, i.e. one can now write x * pf in order to
% scale the @PoleFigure pf by the factor x 
%
% See also
% PoleFigure_index PoleFigure/plus PoleFigure/minus

pf = arg1 .* arg2;
