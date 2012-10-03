function slide=i_demo
% This is a slideshow file for use with playshow.m and makeshow.m
% To see it run, type 'playshow i_demo', 

% Copyright (c) 1984-98 by The MathWorks, Inc.
if nargout<1,
  playshow i_demo
else
  %========== Slide 1 ==========

  slide(1).code={
   '' };
  slide(1).text={
   'Hoi dit is GRINDEMO',...
   'A demo for GRIND'};

  %========== Slide 2 ==========

  slide(2).text={
   'use ZA' };
  slide(2).code={
   'i_use(''za.ini'',1,1,1,0);'};
    %========== Slide 3 ==========

  slide(3).code={
   'null' };
  slide(3).text={
   'null'};
end