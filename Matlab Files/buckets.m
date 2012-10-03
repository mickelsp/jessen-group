%% Buckets Algorithm
%%% Divide independent variable into bins ("buckets"), average dependent
%%% variable values with bins, and determine statistical uncertainty of
%%% averaged value.

%% Revision History
%%% R1 c. 2008, Sarah codes basic algorithm
%%% R2 early 2009, Pascal generalizes algorithm and calculates statistical uncertainty
%%% R3 2009.10.05 Pascal turns algorithm into its own function

%% Common Runtime Errors
%%%1) The thing that causes the most problems is if the the bucket size is
%%%defined poorly. Then, no data ends up in any of the buckets, and the
%%%output vectors are unable to be assigned because they're empty.
%%%2) Another potential problem is if the bucket size is large, but data is 
%%%sparse: the average value gets placed in the center of the bucket, even
%%%though the underlying values are off to the side.  Reduce the bucket
%%%size to alleviate this problem.

%% Function Initialization
%%%INPUTS: vector of independent variable values, vector of dependent
%%%variable values, and the desired bucket size.
%%%OUTPUTS: vector of binned independent variable values, vector of
%%%averaged dependent variable values for each bucket, and statistical
%%%uncertainty of the average dependent variable value.

function [xbuckets,ybuckets,yuncertainty] = buckets(xinput,yinput,bucketsize)
xvector=xinput(:);yvector=yinput(:); %Define input vectors
numlines=size(xvector,1); %How many data points are there?

numbuckets=ceil((max(xvector)-min(xvector))/bucketsize)+1; %Figure out an appropriate number of buckets based on input bucket size; add one to insure inclusion of all points
%%%bucketsize=0.25;   %What it says: the size of the bucket; defined as a function input
bucketoffset=min(xvector); %What minimum value do the buckets begin at?
low=zeros(numbuckets,1);high=zeros(numbuckets,1);counter1=zeros(numbuckets,1);bucketsum=zeros(numbuckets,1); %initialize vectors

%% Create the buckets themselves
for b=1:1:numbuckets %bucket counter
    low(b)=(b-1).*bucketsize+bucketoffset; % establishes lowerbounds
    high(b)=(b).*bucketsize+bucketoffset; %establishes upperbounds
    counter1(b)=0; %establishes counter=0 ---> this variable tells us how many entries there are for a certain time
    bucketsum(b)=0;
end

%% Do the binning (sum all values in the bin)...
for linecount=1:numlines
    for b=1:numbuckets
        if (low(b)<=xvector(linecount)) && (xvector(linecount)<high(b))
            if counter1(b)==0
                bucketsum(b)=yvector(linecount);
            else
                bucketsum(b)=bucketsum(b)+yvector(linecount);
            end
            counter1(b)=counter1(b)+1;
        end
    end
end
%% ...and find the average value in each of the bins
counter2=0; %counts number of bins that actually have values in them
for b=1:numbuckets
    if counter1(b)==0
        bucketsum(b)=0;
    else
        counter2=counter2+1;
        bucketsum(b)=bucketsum(b)/counter1(b); %average number of all data files
        bucketaverage(counter2)=bucketsum(b); %local temporary vector of average of all values in bucket (becomes output y value)
        xaxisbuckets(counter2)=(high(b)+low(b))./2; %local temporary vector of bucket values (becomes output x vector)
    end
end

%% Find the statistical uncertainty for value of each bin
%%%Calculation is for the standard deviation of the mean (Equation 4.14 in John
%%%Taylor's "An Introduction to Error Analysis: The Study of Uncertainties
%%%in Physical Measurements".
numpoints=zeros(counter2,1);sumsquare=zeros(counter2,1);sigmaerror=zeros(counter2,1);
for i=1:counter2 %calculate uncertainty for each bucket that contains values
    counter3=0; %tracks which value we're on within a bucket; reset this count each time we've finished uncertainty for one of the buckets
    for b=1:numbuckets %loop through all possible buckets, filled or unfilled
        if (low(b)<=xaxisbuckets(i)) && (xaxisbuckets(i)<high(b)) %only determine uncertainty when we have a filled bucket
            for linecount=1:numlines %loop through all individual (non-binned) data points
                if (low(b)<=xvector(linecount)) && (xvector(linecount)<high(b)) %only determine uncertainty if individual data point falls in buckets
                    if counter3==0 %when a value is the first in this bucket, start the sum
                        sumsquare(i)=(yvector(linecount)-bucketaverage(i)).^2; %sum of the square of the difference between data point and mean value
                    else %otherwise, add to the existing sum
                        sumsquare(i)=sumsquare(i)+(yvector(linecount)-bucketaverage(i)).^2;
                    end
                    counter3=counter3+1;
                end %end of case structure determining whether unbinned data point falls into the bucket
            end %end of loop through all the individual (non-binned) data points
        end %end of case structure that only determines uncertainty for filled buckets
    end %end of loop through all possible bucket values (b)
    numpoints(i)=counter3;
    sigmaerror(i)=sqrt(1/counter3)*sqrt((1/(counter3-1))*sumsquare(i));
end %of loop through all filled buckets

%% Error handling
%%% It is possible that no data will be found in the specified buckets.
%%% Here, we make sure that buckets.m still outputs something sensible
%%% (empty vectors and an error message) instead of stopping a calling
%%% script from running altogether. 
if counter2==0
    xaxisbuckets=[];bucketaverage=[];sigmaerror=[]; %
    bucketsstatus = 'Bucket number and bucket size are defined such that no data falls into the buckets. Check your input parameter values.'
end

%% Assign values to outputs
xbuckets=xaxisbuckets;
ybuckets=bucketaverage;
yuncertainty=sigmaerror;