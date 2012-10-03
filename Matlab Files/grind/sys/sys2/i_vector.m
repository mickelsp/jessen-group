%i_vector internal function to create a matrix with vectors
function [X, Y, Vect] = i_vector(npoints, iX, Xaxis, iY, Yaxis, N0, relative)
if nargin<7
   relative=0;
end;
global g_grind;
if abs(Xaxis(1)) < 0.00001*abs(Xaxis(2))
   Xaxis(1) = 0.00001*abs(Xaxis(2));
end
if abs(Yaxis(1)) < 0.00001*abs(Yaxis(2))
   Yaxis(1) = 0.00001*abs(Yaxis(2));
end
if (nargin < 6)||isempty(N0)
   N = i_initvar;
else
   N=N0;
end;
findnull=isfield(g_grind,'findnull');
if findnull
   findnull=~isempty(g_grind.findnull.ndx) && ~g_grind.solver.isdiffer;
end;
if findnull
   opt = optimset('fminsearch');
   opt.TolFun = 1E-6;
   opt.TolX = 1E-6;
   opt.MaxFunEval = length(N0) * 8000;
   opt.MaxIter = length(N0) * 8000;
   opt=optimset(opt,'disp','off');
end;
oldX = i_getparvar(iX);
oldY = i_getparvar(iY);
if iX.isext
   oldactX=g_grind.externvars{iX.no}.options.active;
   g_grind.externvars{iX.no}.options.active=0;
end
if iY.isext
   oldactY=g_grind.externvars{iY.no}.options.active;
   g_grind.externvars{iY.no}.options.active=0;
end

try
   nvar = g_grind.statevars.dim;
   Vect = zeros(npoints, npoints, nvar);
   X = zeros(npoints, npoints);
   Y = zeros(npoints, npoints);
   t = 0;
   incrY = (Yaxis(2) - Yaxis(1)) / npoints;
   incrX = (Xaxis(2) - Xaxis(1)) / npoints;
   isdiffer = g_grind.solver.isdiffer;
   for y1 = 1:npoints+1
      py = (y1 - 1) * incrY + Yaxis(1);
      N=i_setparvar(iY,N,py);
      for x1 = 1:npoints+1
         px = (x1 - 1) * incrX + Xaxis(1);
         N=i_setparvar(iX,N,px);
         if findnull
            g_grind.findnull.N0=N;
            Nres = fminsearch('i_findnull', N(g_grind.findnull.ndx), opt);
            g_grind.findnull.N0(g_grind.findnull.ndx)=Nres;
            N=g_grind.findnull.N0;
         end;
         Nres = feval(g_grind.odefile, t + 0.00001, N); %add some small value to avoid
         %recalculation of rednoise
         X(x1, y1) = px;
         Y(x1, y1) = py;
         if isdiffer
            Nres = Nres - N;
         end
         %for v = 1:nvar
         if relative
            Vect(x1, y1, :) = Nres./N;
         else
            Vect(x1, y1, :) = Nres;
         end;
      end;
   end;
   N=i_setparvar(iY,N,oldY);
   i_setparvar(iX,N,oldX);
if iX.isext
   g_grind.externvars{iX.no}.options.active=oldactX;
end
if iY.isext
   g_grind.externvars{iY.no}.options.active=oldactY;
end
   
catch err
 %  err=lasterror;
   N=i_setparvar(iY,N,oldY);
   i_setparvar(iX,N,oldX);
 if iX.isext
   g_grind.externvars{iX.no}.options.active=oldactX;
end
if iY.isext
   g_grind.externvars{iY.no}.options.active=oldactY;
end
  rethrow(err);
end;
