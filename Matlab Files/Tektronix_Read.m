%% Read Data from Tektronix Oscilloscope CSV Files
function [timebase signal numpoints] = Tektronix_Read(filepath) %filepath input should be full file path, not just file name

%% Open File
%filepath = '/Users/work/Desktop/testfile.csv';
fid = fopen(filepath);
rawdata = textscan(fid,'%s %s %f %f %f','delimiter',',');

%% Oscilloscope Settings
numpoints = str2num(rawdata{2}{1}); %[arb] number of data points taken by scope
sample_interval = str2num(rawdata{2}{2}); %[?]
trigger_point = str2num(rawdata{2}{3}); %[s] point at which trigger occurs

vert_scale = str2num(rawdata{2}{9}); %[Volts] vertical scale
vert_offset = str2num(rawdata{2}{10}); %[Volts] how far offset is channel?
horiz_scale = str2num(rawdata{2}{12}); %[s] time scale

%% Extracted Data
timebase = rawdata{4}; %[s] horizontal axis
signal = rawdata{5}; %[volts] signal from channel