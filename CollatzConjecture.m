% Collatz Conjecture

start_val   = 67;

%
V = start_val;
stopFlag = 0;   
k = 1;
while stopFlag == 0
    k = k+1;
    if rem(V(k-1),2) == 0 % even
        V(k) = V(k-1)/2;
    else
        V(k) = V(k-1)*3+1;
    end
    if V(k) == 1
        stopFlag=1;
    end
end
nSteps = k-1;
disp(['# of Steps to Finish: ' num2str(nSteps)]);
plot(V,'lineWidth',1.2)
