function movplot(type, Xs, Ys, Zs, opt)
dim = [];
if nargin == 2
   Ys = [];
   Zs = [];
elseif nargin == 3
   if numel(Ys) == 2
      dim = Ys;
      Ys = [];
   end;
   Zs = [];
end;
plt.type = type;
plt.Xs = Xs;
plt.Ys = Ys;
plt.Zs = Zs;
if isfield(opt,'xlabel')
    plt.xlabel=opt.xlabel;
end;
if isfield(opt,'ylabel')
    plt.ylabel=opt.ylabel;
end;
if isfield(opt,'zlabel')
    plt.zlabel=opt.zlabel;
end;
if ~isempty(Zs)
   % replace 2d functions with 3d variants
   if strcmpi(type,'plot')||strcmpi(type,'line')
      plt.type='plot3';
   elseif strcmpi(type,'bar')
      plt.type='bar3';
   elseif strcmpi(type,'stem')
      plt.type='stem3';
   elseif strcmpi(type,'scatter')
      plt.type='scatter3';
   elseif strcmpi(type,'barh')
      plt.type='bar3h';
   end;
else
   % replace 3d functions with 2d variants
   if strcmpi(type,'plot3')||strcmpi(type,'line3')
      plt.type='plot';
   elseif strcmpi(type,'bar3')
      plt.type='bar';
   elseif strcmpi(type,'stem3')
      plt.type='stem';
   elseif strcmpi(type,'scatter3')
      plt.type='scatter';
   elseif strcmpi(type,'bar3h')
      plt.type='barh';
  end;
end;
plt.dim = dim;
i_movie('movplot', plt);
