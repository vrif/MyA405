function newPress=nudgepress(pressVec)
    newPress=pressVec;
    hit=find(abs(diff(newPress)) < 1.e-8);
    newPress(hit+1)=pressVec(hit) + 1.e-3*pressVec(hit);
end