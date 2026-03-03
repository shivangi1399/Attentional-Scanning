function mergedSac = mergeTwoSaccades(sac1,sac2,maxISI)
%MERGETWOSACCADES(sac1,sac2) Merge Saccade1 and Saccade2
if nargin == 2; maxISI = 0.02; end % in seconds

assert(isa(sac1,'Saccade'),'merge failed, input 1 is not Saccade');
assert(isa(sac2,'Saccade'),'merge failed, input 2 is not Saccade');
assert(~isempty(sac1.parent),'merge failed, input1 has no parent GazeEpoch');
assert(~isempty(sac2.parent),'merge failed, input2 has no parent GazeEpoch');
assert(isequal(sac1.parent, sac2.parent),'merge failed, no common parent GazeEpoch');

parentGaze = sac1.parent;
ind1 = sac1.startIndInParent;
ind2 = sac1.endIndInParent;
ind3 = sac2.startIndInParent;
ind4 = sac2.endIndInParent;

if ind1 <= ind3 % sac1 start before sac2
    isi = parentGaze.time(ind3) - parentGaze.time(ind2);
    assert(isi<maxISI, sprintf('merge failed, isi exceeded maxISI (%f)',maxISI));
    mergedSac = Saccade();
    mergedSac.setParent(parentGaze);
    mergedSac.fetchDataFromParent(ind1,max(ind2,ind4));
    mergedSac.logInfo('merged');
else % sac1 start after sac2
    isi = parentGaze.time(ind1) - parentGaze.time(ind3);
    assert(isi<maxISI, sprintf('merge failed, isi exceeded maxISI (%f)',maxISI));
    mergedSac = Saccade();
    mergedSac.setParent(parentGaze);
    mergedSac.fetchDataFromParent(ind3,max(ind2,ind4));
    mergedSac.logInfo('merged');     
end