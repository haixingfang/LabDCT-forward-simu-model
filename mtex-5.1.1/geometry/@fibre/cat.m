function f = cat(dim,varargin)
% 

% remove emtpy arguments
varargin(cellfun('isempty',varargin)) = [];
f = varargin{1};

fh = cell(size(varargin)); 
fo1 = fh;
fo2 = fh;
for i = 1:numel(varargin)
  f2 = varargin{i};
  if ~isempty(f2)
    fh{i} = f2.h;
    fo1{i} = f2.o1;
    fo2{i} = f2.o2;
  end
end

f.h = cat(dim,fh{:});
f.o1 = cat(dim,fo1{:});
f.o2 = cat(dim,fo2{:});
