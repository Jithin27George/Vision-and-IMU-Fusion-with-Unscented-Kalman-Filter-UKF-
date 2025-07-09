disp(uPrev);
 nDash = 27;
    alpha = 0.001;
    k = 1;
    beta = 2;
    lambdaDash = (alpha^2) * (nDash + k) - nDash;
    disp('lambdaDash   ');
    disp(lambdaDash);
WmDash = [(lambdaDash / (nDash + lambdaDash)), ( ones(1, 2*nDash)*(1 / (2 * (nDash + lambdaDash))) )];
disp(WmDash);
a=ones(1, 2*nDash);
disp(a);

Q = 0.02;
% Calculate Qd
Qd = dt*(eye(12) * Q);
PAug = [covarPrev, zeros(15,12); zeros(12,15), Qd];
disp(size(PAug));
XAug = zeros(nDash, ((2*nDash)+1));
disp(size(XAug));
Xt = zeros(15, ((2*nDash)+1));
disp(size(Xt));

sqrtPAug=chol(PAug);


uAugPrev = [uPrev; zeros(12,1)];
disp('uAugPrev   ');
disp(size(uAugPrev));
for i = 2:1:(nDash+1)
        XAug(:,i) = uAugPrev + sqrt(nDash + lambdaDash) * sqrtPAug(:,i-1);
        XAug(:,(i+nDash)) = uAugPrev - sqrt(nDash + lambdaDash) * sqrtPAug(:,i-1);
end

for i = 2:1:(nDash+1)
    XAugAdd(:,i) = uAugPrev + sqrt(nDash + lambdaDash) * sqrtPAug(:,i-1);
end
for i = 2:1:(nDash+1)
    XAugSub(:,(i+nDash)) = uAugPrev - sqrt(nDash + lambdaDash) * sqrtPAug(:,i-1);
end
disp('XAug   ');
disp(size(XAug));

disp('XAugAdd   ');
disp(XAugAdd);
disp('XAugSub   ');
disp(XAugSub);




 Gi=[(cos(XAug(6,i))*sin(XAug(5,i)))/(cos(XAug(5,i))*cos(XAug(6,i))^2 + cos(XAug(5,i))*sin(XAug(6,i))^2), (sin(XAug(6,i))*sin(XAug(5,i)))/(cos(XAug(5,i))*cos(XAug(6,i))^2 + cos(XAug(5,i))*sin(XAug(6,i))^2), 1;
                                          -sin(XAug(6,i))/(cos(XAug(6,i))^2 + sin(XAug(6,i))^2),                                    cos(XAug(6,i))/(cos(XAug(6,i))^2 + sin(XAug(6,i))^2), 0;
                  cos(XAug(6,i))/(cos(XAug(5,i))*cos(XAug(6,i))^2 + cos(XAug(5,i))*sin(XAug(6,i))^2),              sin(XAug(6,i))/(cos(XAug(5,i))*cos(XAug(6,i))^2 + cos(XAug(5,i))*sin(XAug(6,i))^2), 0];
% ____________________________________________________________________________________________________________________________________________________________
     % R formula (Current body ZYX)
     R = [cos(XAug(6,i))*cos(XAug(5,i)), cos(XAug(6,i))*sin(XAug(4,i))*sin(XAug(5,i)) - cos(XAug(4,i))*sin(XAug(6,i)), sin(XAug(6,i))*sin(XAug(4,i)) + cos(XAug(6,i))*cos(XAug(4,i))*sin(XAug(5,i));
     cos(XAug(5,i))*sin(XAug(6,i)), cos(XAug(6,i))*cos(XAug(4,i)) + sin(XAug(6,i))*sin(XAug(4,i))*sin(XAug(5,i)), cos(XAug(4,i))*sin(XAug(6,i))*sin(XAug(5,i)) - cos(XAug(6,i))*sin(XAug(4,i));
           -sin(XAug(5,i)),                              cos(XAug(5,i))*sin(XAug(4,i)),                              cos(XAug(4,i))*cos(XAug(5,i))];