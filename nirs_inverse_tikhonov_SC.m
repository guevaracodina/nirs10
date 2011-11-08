function out = nirs_inverse_tikhonov_SC(job)
Y = job.Y;
tikh_mask = job.tikh_mask;
alpha = job.alpha;
alpha2 = job.alpha2;
Cp = job.Cp;

switch tikh_mask
    case 'wgmc'
        %masks c1 and c2 (YsegRR)
        m = zeros(size(Y));
        m(Y==1 | Y==2)=1;
        M = m(:);
    case 'image'
        V = spm_vol(Y);
        m = spm_read_vols(V);
        M = m(:);
end

if Cp==0
    out =sparse(diag(alpha*(M+alpha2*(ones(size(M,2),1)-M))));
else
    out =sparse(diag(alpha*Cp*(m(:)+alpha2*(ones(size(M,2),1)-M))));
end
end