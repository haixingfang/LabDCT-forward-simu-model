function varargout = surf(m,cdata,varargin)
%
% Syntax
%
% Input
%
% Output
%
% Options
%
% See also
%

varargin = [m.CS.plotOptions,varargin];
[varargout{1:nargout}] = surf@vector3d(m,cdata,varargin{:},m.CS);
