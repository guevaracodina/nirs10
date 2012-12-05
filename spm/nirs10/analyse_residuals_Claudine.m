function A = analyse_residuals_Claudine%(job)
calcJB = 0;
Su = 1;
s0 = 1;
byCh = 0;
str = 'StatNL';
job.NIRSmat{1} = ['W:\Claudine\SPMDataL\S0' gen_num_str(Su,2) '\' str '\NIRS.mat'];
job.NIRSmatCopyChoice = [];
job.force_redo = 0;
[NIRS newNIRSlocation]=nirs_load(job.NIRSmat{1,1},job.NIRSmatCopyChoice,job.force_redo);
[dir0 fil0] = fileparts(newNIRSlocation);
dir1 = ['W:\Claudine\SPMDataL\ResAnalysis'];
if ~exist(dir1,'dir'), mkdir(dir1); end
fres = fullfile(dir0,['res_Sess' int2str(s0) '.nir']);
NC = NIRS.Cf.H.C.N;
NC = 3*NC/2;
r0 = fopen_NIR(fres,NC);
% for i0=1:NC
% figure; hist(r0(i0,:),100)
% end
if calcJB
    for i0=1:NC
        [h(i0),p(i0),jbstat(i0),critval(i0)] = jbtest(r0(:));
        disp([num2str(h(i0)) '_' num2str(p(i0))])
    end
else
    if byCh
        for c1=1:NC/3
            [ny,nx] = hist(r0(c1,:),100);
            [sigma,mu,A]=mygaussfit(nx,ny,0.2);
            fh0 = figure; hist(r0(c1,:),100); hold on; plot(nx,A*exp(-(nx-mu).^2/(2*sigma^2)),'r');
            legend(['Sigma = ' num2str(sigma)]);
            filen2 = fullfile(dir1,['S0' gen_num_str(Su,2) '_C' gen_num_str(c1,3) '_' str '.tiff']);
            print(fh0, '-dtiffn', filen2);
            close(fh0);
        end
    else
        [ny,nx] = hist(r0(:),100);
        [sigma,mu,A]=mygaussfit(nx,ny,0.2);
        fh0 = figure; hist(r0(:),100); hold on; plot(nx,A*exp(-(nx-mu).^2/(2*sigma^2)),'r');
        legend(['Sigma = ' num2str(sigma)]);
        filen2 = fullfile(dir1,['S0' gen_num_str(Su,2) '_' str '.tiff']);
        print(fh0, '-dtiffn', filen2);
        close(fh0);
    end
end
a=1;

