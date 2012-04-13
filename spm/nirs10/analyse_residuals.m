function A = analyse_residuals%(job)
s_all2 = [1:6 8:11 13 15 17:18 21:22];
Su = s_all2;
AT = 'H'; %Analysis type or flag
s0 = 1; %session number
calcJB = 0;
byCh = 0;
for Idx=Su
    job.NIRSmat{1} = ['D:\Users\Philippe Pouliot\SaidNT\Analysis\S0' gen_num_str(Idx,2) '\Stat' AT '\NIRS.mat'];
    %job.NIRSmat{1} = 'D:\Users\Philippe Pouliot\SaidNT\Analysis\S001\StatJ\NIRS.mat';
    job.NIRSmatCopyChoice = [];
    job.force_redo = 0;
    [NIRS newNIRSlocation]=nirs_load(job.NIRSmat{1,1},job.NIRSmatCopyChoice,job.force_redo);
    [dir0 fil0] = fileparts(newNIRSlocation);
    dir1 = ['D:\Users\Philippe Pouliot\SaidNT\Analysis\ResAnalysis_' AT];
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
        end
    else
        if byCh
            for c1=1:NC/3
                [ny,nx] = hist(r0(c1,:),100);
                [sigma,mu,A]=mygaussfit(nx,ny,0.2);
                fh0 = figure; hist(r0(c1,:),100); hold on; plot(nx,A*exp(-(nx-mu).^2/(2*sigma^2)),'r');
                legend(['Sigma = ' num2str(sigma)]);                
                filen2 = fullfile(dir1,['S0' gen_num_str(Idx,2) '_S' int2str(s0) '_C' gen_num_str(c1,3) '.tiff']);
                print(fh0, '-dtiffn', filen2);
                close(fh0);
            end
        else
            [ny,nx] = hist(r0(:),100);
            [sigma,mu,A]=mygaussfit(nx,ny,0.2);
            fh0 = figure; hist(r0(:),100); hold on; plot(nx,A*exp(-(nx-mu).^2/(2*sigma^2)),'r');
            legend(['Sigma = ' num2str(sigma)]);
            
            filen2 = fullfile(dir1,['S0' gen_num_str(Idx,2) '_S' int2str(s0) '.tiff']);
            print(fh0, '-dtiffn', filen2);
            close(fh0);
        end
    end
    a=1;
end
