function xBF = nirs_specify_bases(job)
if  job.units == 0
    xBF.UNITS = 'scans';
elseif  job.units == 1
    xBF.UNITS = 'secs';
end
str_hrf = 'hrf';
str_eHRF = 'eHRF';
%specify the basis functions of the HRF to be used
str0 = '';
if strcmp(fieldnames(job.bases),str_hrf)
    str0 = str_hrf;
else
    if strcmp(fieldnames(job.bases),str_eHRF) 
        str0 = str_eHRF;
    end
end

if ~isempty(str0)
    xBF.name0 = str0;
    if all(job.derivs == [0 0])
        xBF.name = str0;
    elseif all(job.derivs == [1 0])
        xBF.name = [str0 ' (with time derivative)'];
    elseif all(job.derivs == [1 1])
        xBF.name = [str0 ' (with time and dispersion derivatives)'];
    else
        disp('Unrecognized hrf derivative choices.')
    end
else
    nambase = fieldnames(job.bases);
    if ischar(nambase)
        nam=nambase;
    else
        nam=nambase{1};
    end
    switch nam,
        case 'fourier',
            xBF.name = 'Fourier set';
        case 'fourier_han',
            xBF.name = 'Fourier set (Hanning)';
        case 'gamma',
            xBF.name = 'Gamma functions';
        case 'fir',
            xBF.name = 'Finite Impulse Response';
        otherwise
            error('Unrecognized hrf derivative choices.')
    end
    xBF.length = job.bases.(nam).length;
    xBF.order  = job.bases.(nam).order;
end