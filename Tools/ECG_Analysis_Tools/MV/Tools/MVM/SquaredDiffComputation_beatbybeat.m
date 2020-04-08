function [NTWDseries TWDseries] = SquaredDiffComputation_beatbybeat(Complexes)
% [NTWDseries TWDseries] = SquaredDiffComputation_beatbybeat(Complexes)
%   OVERVIEW:   Computed squared difference between consecutive QRS
%               complexes with and without dynamic time warping.
%
%   INPUT:      MANDATORY:
%               Complexes       : 2D array of dimension number of complexes by median length of complexes
%
%               
%                                                
%                
%   OUTPUT:     
%               NTWDseries      : Squared difference between consecutive QRS complexes without time warping
%
%               TWDseries       : Squared difference between consecutive QRS complexes with time warping                            
%
%	REPO:       
%       https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox
%   ORIGINAL SOURCE AND AUTHORS:     
%       Written by Ismail Sadiq    
%	COPYRIGHT (C) 2019
%   LICENSE:    
%       This software is offered freely and without warranty under 
%       the GNU (v3 or later) public license. See license file for
%       more information. The license may be found in
%       the Documents folder of the Physionet-Cardiovascular-Signal-Toolbox.
%

noofconsecutivecomplexes = size(Complexes,1);

% time warped squared diff
TWDseries(1:noofconsecutivecomplexes) = 0;
% non time warped squared diff
NTWDseries(1:noofconsecutivecomplexes) = 0;

for index = 1:noofconsecutivecomplexes-1
    
    complex1 = Complexes(index,:);    
    complex2 = Complexes(index+1,:);    
    
    % compute squared difference between each complex and median complex w/o
    % dtw
    SD = sum((complex1 - complex2).^2);
    SD = SD ./ length(complex1);
    NTWDseries(index) = SD;
    
    % compute squared difference between each complex and median complex with dtw
    [dist,ix,iy] = dtw(complex1, complex2);
    templt1 = complex1(ix);
    templt2 = complex2(iy);
    
    MD = sum((templt1 - templt2).^2);
    MD = MD ./ length(templt1);
    
    TWDseries(index) = MD;
    
end

end

