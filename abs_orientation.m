function [s R t] = abs_orientation(x,y)
% y = sR(x) + t
% find s,  R,  t

N = size(x,2);

% centroid
cntrX = (sum(x')/N)';
cntrY = (sum(y')/N)';

% new coordinate
xNew = x - cntrX*ones(1,N);
yNew = y - cntrY*ones(1,N);

% make 'NN' matrix
NN = zeros(4,4);
for(iterData = 1:N)
    NN = NN + quaternion2(xNew(:,iterData))'*quaternion1(yNew(:,iterData));
end;

% find maximun eigenvector
[eigv eigd] = eig(NN);
qtn = eigv(:,4);

% estimate R
R = qtn2rtm(qtn);

% find scale factor
s = 0;
for iterData = 1:N
    ry_dot_R_rx = sum(yNew(:,iterData) .* (R * xNew(:,iterData)));
    s = s + ry_dot_R_rx;
end;
s = s / sum(xNew(:).^2);

% find translation factor
t = cntrY - s * R*cntrX;