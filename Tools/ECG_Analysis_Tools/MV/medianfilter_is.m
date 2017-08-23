function [ signal_filter2 ] = medianfilter_is( signal, fs )
%medianfilter_is: removes baseline wander from signal using median filter
%   Detailed explanation goes here

orderfilter1 = floor(0.2 * fs);
orderfilter2 = floor(0.6 * fs);

signal_filter1 = medfilt1(signal, orderfilter1);
signal_filter2 = medfilt1(signal_filter1, orderfilter2);

end