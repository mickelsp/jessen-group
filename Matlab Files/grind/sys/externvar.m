%EXTERNVAR   Get the value of an external variable
%   This function can evaluate an external variable. Use defextern to 
%   define external variables. The function interpolates the external (file) data. 
%   The used method is linear interpolation. 
%
%
%   Usage:
%   EXTERNVAR(varnr,default,t) - the function has three parameters: (1) variable number (see
%   g_grind.externalvars. 2) default, a value which is used outsite the range of the data
%   (3) t, current time or a vector with time.
%   EXTERNVAR(varnr,t) - if no default value is given, the function returns NaN outside
%   the range of the data.
%   EXTERNVAR(name,t) - You can also use the name of the variable.
%
%
%   See also defextern, setdata, loaddata, model, rednoise

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:26 $
function V = externvar(varnr, default, at)
global g_grind g_t;
if nargin == 2
   at = default;
   default = NaN;
elseif (nargin==1)&& ischar(varnr)
   i=1;
   while i<=length(g_grind.externvars)&&~strcmp(varnr,g_grind.externvars{i}.name)
      i=i+1;
   end;
   varnr=i;
   if i<=length(g_grind.externvars)
      V=externvar(varnr, g_grind.externvars{varnr}.default,g_t);
      return;
   else
      error('GRIND:externvar:UnknownName','Unknown name');
   end;
elseif nargin==0
   error('GRIND:externvar:TooFewArgs','EXTERNVAR: too few arguments of the function');
end;
options=g_grind.externvars{varnr}.options;
data=g_grind.externvars{varnr}.data;
if ~options.active
    V= evalin('base',g_grind.externvars{varnr}.name);
   return;
end;
s = size(data);
if ((s(1) == 2) || (s(1) == 1)) && (s(2) > 2)
   data = data';
   s = size(data);
end;
if options.tofloor
   at=floor(at);
end;
if s(2) == 1
   if options.cycle
      at = mod(at, s(1));
   end;
   if length(at) > 1
      V = interp1((1:s(1))', data, at);
   else
      if (at + 1 > s(1)) || (at < 0)
         V = NaN;
      else
         remain = rem(at + 1, 1);
         V = data(floor(at + 1)) * (1 - remain) + data(ceil(at + 1)) * remain;
      end;
   end;
elseif s(2) == 2
   if options.cycle
      t1=data(1,1); tend=data(end,1);
      at = mod(at, tend - t1)+t1;
   end;
   if isempty(at)
      V=[];
   elseif length(at) > 1
      V = interp1(data(:, 1), data(:, 2), at);
   else
      V = parlookup(data, at,0,0);
   end;
elseif s(1) == 0
   V = ones(size(at)) * NaN;
else
   error('GRIND:externvar:dataerror','EXTERNVAR: data matrix not correct')
end;
V(isnan(V)) = default;
