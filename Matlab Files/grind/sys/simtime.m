%SIMTIME   Set the simulation duration
%   Enter the following information:
%   - Start time - sets the start time (can also be assigned as t =..)
%   - Number of steps to simulate - number of time steps for default runs
%   - Number of steps for output (leave empty for maximal) - you can reduce
%     the number of steps in the output, for example as Poincare sections for
%     yearly cycles. Leave this value empty for maximal output.
%
%    Usage:
%    SIMTIME - enter the data interactively.
%    SIMTIME T STEPS  - sets the starting time to T and the number of time steps to STEPS.
%    SIMTIME T STEPS SOUT - sets the number of steps for output to SOUT also.
%
%
%   See also ru, time

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:29 $
function simtime(at, ndays, tstep)
global t g_grind;
if nargin >= 1
   t = i_checkstr(at);
end
if nargin >= 2
   g_grind.ndays = i_checkstr(ndays);
end;
if nargin == 3
   g_grind.tstep = i_checkstr(tstep);
end;
if nargin == 0
   i_parcheck;
   prompt = {'Start time:', ...
      'Number of time steps to simulate', ...
      'Number of values for output (leave empty for maximal)'};
   answer = cell(3);
   answer{1} = num2str(t);
   answer{2} = num2str(g_grind.ndays);
   if ~isnan(g_grind.tstep)
      answer{3} = num2str(g_grind.tstep);
   else
      answer{3} = '';
   end;
   answer = inputdlg(prompt, 'Simulator time', 1, answer);
   if ~isempty(answer)
      t = i_checkstr(answer{1});
      g_grind.ndays = i_checkstr(answer{2});
      if ~isempty(strtrim(answer{3}))
         g_grind.tstep = i_checkstr(answer{3});
      else
         g_grind.tstep = NaN;
      end;
   end;
end
