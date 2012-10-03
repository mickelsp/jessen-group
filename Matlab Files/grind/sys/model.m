%MODEL   Open/create/edit a model
%   Use the upper panel in the window to enter differential 
%   equation or difference equations (or iteration maps).
%   In the lower panel default values of the parameters and initial
%   values of the state variables are entered.  
%   Furthermore default commands can be entered here (for instance
%   axis limits, see ax).
%   There are no restrictions for the number of parameters, functions and equations.
%
%   Note that parameters and state variables are case-sensitive.
%   Parameter names can be all alphanumeric names except the
%   following reserved words:
%     t   - initial time.
%     pi  - pi=3.1416
%     inf - Infinity
%     Inf - Infinity
%     nan - NaN (not-a-number)
%     NaN - NaN (not-a-number)
%     eps - Floating point relative accuracy eps=2.2204e-016 
%
%
%   Examples:
%   Model equations:
%   cons=Z*A/(h+A)  - "function" with temporal results that can be used in the 
%   other equations (see funcs) Such function cannot take arguments.
%   N(t + 1) = N(t) * r * (1 - N(t) / K)  - Use (t+1) and (t) for difference equations 
%   Higher order equations (N(t + 2) ) are not allowed, but can be written as first
%   order equations.
%   N' = N * r * (1 - N / K) - Use ' for differential equations with respect of time (dN/dt).
%   X(1:4)'=X .* (gamma - V * X) - You may also use such matrix notation. The length 
%   of the vector of state variables must be set by use of a colon. *. is a product of 
%   arrays, * is a matrix product. See MATLAB manuals.
%   X(1:10,1:10)(t)= f(X(t-1)) - State variables can also be matrices (in this example a 10x10 
%   matrix).  There are special commands for viewing the matrices (viewcells).
%   function res=monod(x,h);
%   res=x/(x+h);
%   return; - More complex functions with arguments. All code between function 
%   and return is copied into a function file. You may use any MATLAB statement here,
%   but matrix multiplication is not supported. Finish your function always with "return".
%
% 
%   Default values of the parameters/Commands:
%   N=0.01;  - assign initial values and default parameters.
%   r=0.5;
%   ax x N [0 10];
%   (the semicolon is not required, but suppresses unnecessary
%   output)
%   setevent('simpleevent',0,'A=A+1',30); 
%  
%
%   Usage:
%   MODEL - Opens with the current model
%   MODEL AMODEL - Opens AMODEL
%   MODEL -clear - Clears the current model from memory (and deletes currently temporary files)
%
%   See also use, savemodel, vismod, lag, rednoise, modelpanel

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:27 $
function model(varargin)
global g_grind;
%if the model name has spaces it can be read as separate arguments
if (nargin==1) && strncmpi(varargin{1}, '-c', 2)
   finishgrind;
   if ~isempty(g_grind)&&isfield(g_grind,'statevars');
      if g_grind.statevars.vector
         s = sprintf('%s ', g_grind.pars{:}, g_grind.statevars.vectnames{:});
      else
         s = sprintf('%s ', g_grind.pars{:}, g_grind.statevars.names{:});
      end;
      s = sprintf('clear global %s t g_*', s);
      evalin('base', s);
      rmpath([grindpath filesep 'sys2']);
      disp('Cleared model from memory and deleted temporary files');
   else
      if exist('i_use', 'file')==2
         rmpath([grindpath filesep 'sys2']);
      end;
      warning('GRIND:model:noopenfile','No GRIND model is open');
   end;
   return;
end;
if ~exist('i_use', 'file')
  addpath([grindpath filesep 'sys2']);
end;
amodel = strtrim(sprintf('%s ', varargin{:}));
if isempty(amodel)&&~isempty(g_grind)&&(length(g_grind.scheme)>10)
    warning('GRIND:model:modelhasscheme','The current model has a Forrester scheme, and should be edited with <a href="matlab:vismod">vismod</a>');
   % vismod;
   % return;
end;
i_use(amodel, 1, 0);


