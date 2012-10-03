%SETDIMENSION   Change the dimension of a (state) variable
%   This function changes the dimension of a vector/matrix state variable or parameter. It changes
%   state variables dynamically but it doesn not check for neccessary changes in vector
%   parameters. Usually this function is used in scripts/programs.
%   Note that if this is done during a run, model results can get lost if the dimension is shrinking. 
%
%
% Usage:
%    SETDIMENSION('VAR',NROWS,NCOLS) - set the dimension of variable/parameter VAR to NROWS and NCOLS 
%    SETDIMENSION('VAR',NROWS) - assumes NCOLS=1, or NCOLS=NROWS if the variable is already square 
%
%  See also model

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:28 $

% SETDIMENSION - change the dimension of a state variable
%   This function changes the dimension of a vector/matrix state variable, however it doesn not check for neccessary
%   changes in vector parameters. Usually this function is used in scripts/programs.
%   Note that if this is done during a run, data can get lost if the dimension is shrinking. Output for time plots
%   is set to the default values;
%
% Usage:
%    SETDIMENSION('VAR',NROWS,NCOLS) - set the dimension of variable VAR to NROWS and NCOLS
%    SETDIMENSION('VAR',NROWS) - assumes NCOLS=1
%
function res=setdimension(avar, dim1, dim2)
siz=evalin('base',sprintf('size(%s)',avar));
if nargin < 3
   if siz(1)==siz(2)
      dim2=dim1;
   else
      if siz(1)~=1
         dim2 = 1;
      else
         dim2=dim1;
         dim1=1;
      end;
   end;
end
if nargin < 2
   error('GRIND:setdimension:ArgError','Not enough arguments for setdimension');
end;
dim1 = i_checkstr(dim1);
dim2 = i_checkstr(dim2);
global g_Y g_grind g_func;
if isempty(g_grind)||~g_grind.statevars.vector
   setdimvar(avar,dim1,dim2,siz);
   if nargout>0
      res=evalin('base',avar);
   end;
   return;
end;
ivect = 0;
for i = 1:size(g_grind.statevars.vectnames, 2)
   if strcmp(avar, g_grind.statevars.vectnames{i});
      ivect = i;
      break;
   end;
end;
if ivect == 0
   setdimvar(avar,dim1,dim2,siz);
    if nargout>0
      res=evalin('base',avar);
   end;
   return;
end;
gfunc = g_func;
ggrind = g_grind;
try
   %erase the function results as they might get another dimension
   g_func = [];
   nams = g_grind.funcnames.names;
   g_grind.funcnames = [];
   g_grind.funcnames.names = nams;
   diff = dim1 * dim2 - g_grind.statevars.dims{ivect}.dim1 * g_grind.statevars.dims{ivect}.dim2;
   g_grind.statevars.dims{ivect}.dim1 = dim1;
   g_grind.statevars.dims{ivect}.dim2 = dim2;
   g_grind.statevars.dims{ivect}.to = g_grind.statevars.dims{ivect}.to + diff;
   g_grind.statevars.dim = g_grind.statevars.dim + diff;
   for i = ivect + 1:size(g_grind.statevars.vectnames, 2)
      g_grind.statevars.dims{i}.from = g_grind.statevars.dims{i}.from + diff;
      g_grind.statevars.dims{i}.to = g_grind.statevars.dims{i}.to + diff;
   end;
   g_grind.statevars.names = {};
%  d = 0;
%   for j = 1:length(g_grind.statevars.vectnames)
%      if g_grind.statevars.dims{j}.dim2 == 1
%         for i = 1:g_grind.statevars.dims{j}.dim1
%            g_grind.statevars.names{d + i} = sprintf('%s(%d)', g_grind.statevars.vectnames{j},  i);
%         end;
%         d = d + g_grind.statevars.dims{j}.dim1;
%      else
%         for k  = 1:g_grind.statevars.dims{j}.dim1
%            nbase = d + (k - 1) * g_grind.statevars.dims{j}.dim1;
%            for i = 1:g_grind.statevars.dims{j}.dim1
%               g_grind.statevars.names{nbase + i} = sprintf('%s(%d,%d)',g_grind.statevars.vect2names{j},i,k);
%            end;
%         end;
%         d = d + g_grind.statevars.dims{j}.dim1 * g_grind.statevars.dims{j}.dim2;
%      end;
%   end;
   if ~isempty(g_Y)  && (size(g_Y, 2) < g_grind.statevars.dim)
      g_Y(size(g_Y, 1), g_grind.statevars.dim) = 0;
   end;
   setdimvar(avar,dim1,dim2,siz);
   if nargout>0
      res=evalin('base',avar);
   end;
catch err
%   err=lasterror;
   g_grind = ggrind; %#ok
   g_func = gfunc;  %#ok
   rethrow(err);
end;

%out('-defaults');
function setdimvar(avar,dim1,dim2,siz)
   diff=dim1*dim2-prod(siz);
   if diff < 0
      if (dim1<=siz(1))&&(dim2<=siz(2))
         evalin('base',sprintf('%s=%s(1:%d,1:%d);',avar,avar,dim1,dim2));
      else
         evalin('base',sprintf('%s=reshape(%s(1:%d),%d,%d);',avar,avar,dim1*dim2,dim1,dim2));
      end;
   elseif diff == 0
      evalin('base',sprintf('%s=reshape(%s,%d,%d);',avar,avar,dim1,dim2));
   elseif (dim1<=siz(1))&&(dim2<=siz(2))
      evalin('base',sprintf('h=%s;%s=zeros(%d,%d)+mean(mean(%s));%s(1:%d,1:%d)=h;',avar,avar,dim1,dim2,avar,avar,siz(1),siz(2)));    
   else
      evalin('base',sprintf('h=%s;%s=zeros(%d,%d)+mean(mean(%s));%s(1:%d)=h;',avar,avar,dim1,dim2,avar,avar,siz(1)*siz(2)));    
  end;



