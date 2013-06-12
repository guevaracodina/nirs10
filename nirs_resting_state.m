function nirs_resting_state
%Pil04_ED_restingState
load('Pil04_ED_restingState.nirs','-mat');
figure; plot(t,d)
%Bandpass filter data from 0.008 to 0.08 Hz
fs = 1/(t(2)-t(1));
fl = 0.008; fh = 0.08; order = 4;
fls = 0.004; fhs = 0.12; Rp = 3; Rs = 20; %stopband frequencies and attenuations
%frequencies normalized to Nyquist
wn = [2*fl/fs 2*fh/fs];
ws = [2*fls/fs 2*fhs/fs];
[nl wnl] = buttord(wn(1),ws(1),Rp,Rs);
[nh wnh] = buttord(wn(2),ws(2),Rp,Rs);
[nbp wnbp] = buttord(wn,ws,Rp,Rs);
[nbp1 wnbp1] = cheb1ord(wn,ws,Rp,Rs);
[nbp2 wnbp2] = cheb2ord(wn,ws,Rp,Rs);
[nbpe wnbpe] = ellipord(wn,ws,Rp,Rs);
[z,p,k] = butter(nbp,wnbp,'bandpass');
r_pp = 3;
[z1,p1,k1] = cheby1(nbp1,r_pp,wnbp1,'bandpass');

[bb,ab] = butter(nbp,wnbp,'bandpass');
db = filtfilt(bb,ab,d);
[bc1,ac1] = cheby1(nbp1,r_pp,wnbp1,'bandpass');
dc1 = filtfilt(bc1,ac1,d);

%figure; ds = 1; plot(d(:,ds),'r'); hold on; plot(dc1(:,ds),'k');
figure; plot(t,dc1);
cor = corr(dc1);
figure; imagesc(cor)
%l2 = 850 nm
l2 = sort([5:8:128 6:8:128 7:8:128 8:8:128]);
dcl2 = dc1(:,l2);
corl2 = corr(dcl2);
figure; plot(t,dcl2);
figure; imagesc(corl2)
anticor = dcl2(:,[3 4 20 24 28]);
selc = l2([3 4 20 24 28])
ml(selc,:)
figure; plot(t,anticor);