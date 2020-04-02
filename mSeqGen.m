function m = mSeqGen
G = [8 7 6 1 0];
mG = comm.PNSequence('Polynomial', G, 'InitialConditions', 1, ...
    'SamplesPerFrame', 2^8-1, 'Mask', de2bi(1,8));
mSeq = mG();

%M Sequence using Fibonacci Implementation
mLength = 255;
m = zeros(1,mLength); % initialize m sequence vector
s = [1 0 0 0 0 0 0 0]; % intial state
for i = 1:mLength
    x =mod(s(1)+s(6)+s(7)+s(8),2);
    s = circshift(s,1);
    s(1) = x;
    m(i) = s(8);
end

m1Flip = flip(m);
for i = 1:length(m)
    check(i) = sum(m1Flip == transpose(mSeq));
    if check(i) == 255
        break
    end
    m1Flip = circshift(m1Flip,1);
end
verify = max(check) == 255;
m1Flip == mSeq';

m=m1Flip';

end