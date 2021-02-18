%This function is a wrapper for the basic EOF library function to add a
%couple of extra features - particularly, guaranteeing a descending sort
%for the EOFs and providing normalised as well as raw eigenvalues.
function [Eigenvals,EigenvalsNormd,EOFs,ECs,TruncnErr] = EOFWrapper(Data,numEOFs)

[Eigenvals,EOFs,ECs,TruncnErr] = EOF(Data,numEOFs);
[Eigenvals,eigenSortList] = sort(Eigenvals,'descend');
EOFs = EOFs(:,eigenSortList);
EigenvalsNormd = Eigenvals/sum(Eigenvals);

end