function odf = uminus(odf)
% superposeing two ODFs
%
% overload the - operator, i.e. one can now write - @ODF
%
% See also
% ODF_index ODF/mtimes

odf.weights = -odf.weights;

