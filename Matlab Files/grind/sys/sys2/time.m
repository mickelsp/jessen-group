%TIME   Create a time series plot
%   Show the data of the last run in a time plot. On the y-axis are the
%   variables or functions that are determined by OUT. It is possible to update
%   several figures simulaneously, see out how to define this. If any parameters or
%   initial conditions have been changed, first the model is run for g_grind.ndays.
%
%   Usage:
%   TIME - runs the model for g_grind.ndays time steps.
%   TIME N - (or time(N)) runs the model for N time steps.
%   TIME -OPT runs the model using the option OPT. Valid options are (short version 
%          between brackets):
%         -add (-a) continue with the previous run.
%         -hold (-h) do not overwrite previous plot.
%         -file=afile (-f)  write the data to a tab delimited file (if the filename is 
%          omitted, a saveas dialog appears.
%         -nocheck (-n)  do not check for changed option, just show last results.
%         -paranal (-p1)  show the results of the last run of PARANAL.
%         -paranal2 (-p2)  show the results of the last two runs of PARANAL
%         -out (-o) select output for current time plot.
%         -run (-r) rerun the model always.
%         -silent (-s) do not plot the results.
%
%   TIME  N -OPT1 -OPT2 - combine options and time steps 
%   RES = time('N','-OPT') - save the results to the matrix RES 
%    
%
%   See also ru, out, simtime, addmode, outfun

%   Copyright 2012 WUR
%   Revision: 1.1.8 $ $Date: 15-Mar-2012 10:05:29 $
function [outmat, leg] = time(varargin)
if nargout>0
   [outmat, leg] = i_time(varargin{:});
else
   i_time(varargin{:});
end;
