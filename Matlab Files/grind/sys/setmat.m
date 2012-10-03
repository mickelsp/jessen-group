%SETMAT   Set the values of a matrix (disturbance/gradient)
%   Use this command to set easily the values of a matrix.
%   
%   Usage:
%   SETMAT - opens a window to set the values
%   SETMAT A VAL - sets all values of the matrix A to VAL
%   SETMAT A [L1,L2] SIZE VAL - sets in matrix A at [L1,L2] a square of SIZE cells to VAL
%   SETMAT A -centre SIZE VAL - disturbance VAL,SIZE at the centre
%   SETMAT A -random SIZE VAL - disturbance VAL,SIZE at a random position
%   SETMAT A [L1,L2] [S1,S2] VAL - sets a rectancle of S1xS2 to VAL
%   SETMAT A LOC SIZE [V1,V2] - sets outside the disturbance to V1, inside to V2
%   SETMAT A -xgradient [V1,V2] - sets an x gradient from V1 to V2
%   SETMAT A -ygradient [V1,V2] - sets an y gradient from V1 to V2
%   SETMAT A -uniform [V1,V2] - draws at random values between V1 and V2 (default [0 1])
%   
%

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:29 $
function [A, locs] = setmat(var, loc, siz, value, num)
global g_grind;
if (nargin == 2) && (~ischar(loc) || (loc(1) ~= '-'))
   value = loc;
   if ~ischar(value)
      value = num2str(value);
   end;
   loc = '-set';
   siz=9;
else
   if nargin < 5
      num = 1;
   else
      num = i_checkstr(num);
   end;
   if nargin < 3
      siz = 9;
   else
      siz = i_checkstr(siz);
   end;
   if nargin < 4
      value = '0';
   else
      if ~ischar(value)
         value = num2str(value);
      end;
   end;
   if nargin < 2
      loc = '-centre';
   else
      if ~(ischar(loc) && (loc(1) == '-'))
         loc = i_checkstr(loc);
      end;
   end;
   if nargin < 1
      if ~isempty(g_grind)
         var = g_grind.statevars.vectnames{1};
      else
         var = '';
      end;
      if ~isempty(g_grind) && isfield(g_grind, 'setmat')
         answer = g_grind.setmat;
      else
         answer = {var, loc, num2str(siz), value, num2str(num)};
      end;
      prompt={'Name of matrix variable','Location of disturbance/option ([x y] -centre -random)',...
         'Size of disturbance (number of cells) or [nx ny]','New value of disturbance or [outside inside]', 'Divide area in n parts (only with -random)'};
      answer = inputdlg(prompt, 'set matrix', 1, answer);
      if isempty(answer)
         disp('Cancelled');
         return;
      end;
      var = answer{1};
      loc = answer{2};
      if (loc(1) ~= '-')
         loc = i_checkstr(loc);
      end;
      siz = i_checkstr(answer{3});
      value = answer{4};
      num = i_checkstr(answer{5});
      if ~isempty(g_grind)
         g_grind.setmat = answer;
      end;
   end;
end;
siz=round(siz);
if ischar(var)
   E = str2num(var); %#ok
   if isempty(E)
      E = evalin('base', var);
   end;
   n = size(E);
else
   E = var;
   var = '';
   n = size(E);
end;
if ischar(loc)
   if strncmpi(loc, '-c', 2)
      loc = round([n(1), n(2)] ./ 2);
      E = setdisturbance(E, n, loc, siz, value);
   elseif strncmpi(loc, '-xh', 3)
      siz=size(E);
      E(:,1:round(siz(1)/2))=i_checkstr(value);
    elseif strncmpi(loc, '-yh', 3)
      siz=size(E);
      E(1:round(siz(1)/2),:)=i_checkstr(value);
    elseif strncmpi(loc, '-r', 2)
      sizes = intdivide(siz, num);
      found=0;
      % set num non-overlapping not bordering random locations 
      Etmp=E;
      Eopt=[];locsopt=[];
      maxsize=0;
      j=0;
      while ~found
         nuls=zeros(size(E));
         locs=[0 0];
         E=Etmp;
         for i = 1:length(sizes)
            loc = floor([rand(1) * n(1), rand(1) * n(2)]) + 1;
            locs(i,:)=loc;
            E = setdisturbance(E, n, loc, sizes(i), value);
            nuls = setdisturbance(nuls, n, loc, sizes(i), '1');
         end;
         s=powerlaw(nuls);
         ssum=sum(sum(nuls));
         if ssum>maxsize
            Eopt=E;
            locsopt=locs;
            maxsize=ssum;
         end;
         j=j+1;
         found=((ssum==siz) && (max(s)<=max(sizes))) || (j>200);
     end;
     if j>200
        warning('GRIND:setmat:distoverlap','Could not create non-overlapping "disturbances"');
        E=Eopt;
        locs=locsopt;
     end;
   elseif strncmpi(loc, '-u', 2)||strncmpi(loc, '-n', 2)
      if isempty(siz)
         siz = i_checkstr(value);
      end;
      if length(siz) == 1
         if nargin == 3
            siz = [0 siz];
         else
            siz = [0 1];
         end;
      end;
      if isnan(siz(1)) %mean or minimum is NaN  - > add noise to current value
         m = E;
      else
         m = siz(1);
      end;
      if strncmpi(loc, '-u', 2)
         E = rand(size(E)) * siz(2) + m;
      else
         E = randn(size(E))  * siz(2) + m;
         E(E < 0) = 0;
      end;
   elseif strncmpi(loc, '-x', 2)||strncmpi(loc, '-y', 2)
      if isempty(siz)
         siz = i_checkstr(value);
      end;
      if length(siz) == 1
         siz = [0 1];
      end;
      [X,Y] = meshgrid(siz(1):(siz(2) - siz(1)) / (n(1) - 1):siz(2), siz(1):(siz(2) - siz(1)) / (n(2) - 1):siz(2));
      if strncmpi(loc, '-x', 2)
         E = X;%assignin('base', var, X); plotmat(X,var);
      else
         E = Y;%assignin('base', var, Y); plotmat(Y,var);
      end;
   elseif strncmpi(loc, '-s', 2)
      E = ones(size(E)) * evalin('base', value);
      %      assignin('base', var, E); plotmat(E,var);return;
   end;
else
   E = setdisturbance(E, n, loc, siz, value);
end;

if ~isempty(var)
   assignin('base', var, E);
end;
if nargout > 0
   A = E;
else
   plotmat(E, var);
end;

function E = setdisturbance(E, n, loc, siz, value)
if siz > 0
   [xx, yy] = getdistloc(siz);
   xx = xx + loc(1);
   xx(xx <= 0)=xx(xx <= 0) + n(1);
   xx(xx > n(1)) = xx(xx > n(1)) - n(1);
   yy = yy + loc(2);
   yy(yy <= 0)=yy(yy <= 0) + n(2);
   yy(yy > n(2)) = yy(yy > n(2)) - n(2);
   val = evalin('base', value);
   if length(val) == 2
      E = zeros(size(E)) + val(1);
      val = val(2);
   end;
   E(sub2ind(size(E), xx, yy)) = val;
end;

function plotmat(E, var)
global g_grind;
i_makefig('setmat');
if min(size(E)) == 1
   plot(1:length(E), E);
   if ~isempty(g_grind)
      title(sprintf('initial conditions %s', var));
      i_plotdefaults;
   else
      title(sprintf('matrix %s', var));
   end;
   xlabel('Row nr');
else
   pcolor(E);
   if ~isempty(g_grind)
      colormap('i_grindcolormap');
      i_plotdefaults;
      title(sprintf('initial conditions %s', var));
   else
      title(sprintf('matrix %s', var));
   end;
   daspect([1 1 1]);
   set(gca,'drawmode','fast','ydir','reverse');
   xlabel('Column');
   ylabel('Row');
   shading flat;
   colorbar;
   i_plotdefaults;
end;

%get x's and y's of a disturbance of n cells (round zero)
%keep square as much as possible
function [Xs, Ys] = getdistloc(n)
if length(n) == 2
   if sum(abs(n - round(n))) > 0.05
      [Xs, Ys] = getdistloc(prod(n));
   else
      [Xs, Ys] = meshgrid(-floor(n(1) / 2):n(1)-floor(n(1) / 2)-1,...
         -floor(n(2) / 2):n(2)-floor(n(2) / 2)-1);
      Xs = Xs(:)';
      Ys = Ys(:)';
   end;
else
   n = round(n);
   root = floor((floor(sqrt(n)) + 1) / 2) * 2 - 1;
   siz = (root - 1) / 2;
   [Xs, Ys] = meshgrid(-siz:siz, -siz:siz);
   borderX = -[ones(1, siz * 2 + 1) * (-siz-1) (-siz-1:1:siz + 1)  ones(1, siz * 2 + 1) * (siz + 1) (siz + 1:-1:-siz) (-siz-1)];
   borderY = [siz:-1:-siz-1 ones(1, siz * 2 + 1) * (-siz-1) (-siz-1:1:siz + 1) ones(1, siz * 2 + 2) * (siz + 1)];
   Xs = [Xs(:)' borderX(1:n - length(Xs(:)))];
   Ys = [Ys(:)' borderY(1:n - length(Ys(:)))];
end;

function A = intdivide(n, parts)
h = 0;
part = floor(n / parts);
addh = n / parts - part;
A = zeros(1,parts);
for i = 1:parts
   h = h + addh;
   if h > 0.99
      h1 = floor(h + 0.001);
      h = h - h1;
   else
      h1 = 0;
   end;
   A(i) = part + h1;
end;

%function loc=stata(i,num)
%if num<4
%   loc=rand(1,2);
%elseif num<8 
%   loc(2)=rand(1);
%   if i<num/2
%      loc(1)=rand(1)*0.5;
%   else
%      loc(2)=rand(1)*0.5+0.5;
%   end;
%elseif num<12
%   loc=rand(1,2);
%   if i>num/2
%      loc(1)=loc(1)+0.5;
%   end;
%   if (i>num/4)&(i<3*num/4)
%      loc(2)=loc(2)+0.5;
%   end;
%end;



