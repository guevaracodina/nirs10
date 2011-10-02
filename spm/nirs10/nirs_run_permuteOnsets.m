function out = nirs_run_permuteOnsets(job)
%Permute onsets if required
doPermute = 0;
fOnsetOK = 0;
if ~isempty(job.onset_files)
    nf = size(job.onset_files,1);
    if nf < 2
        disp('Nothing to permute, exiting');
        out = [];
    else
        if length(job.NIRSmat) == 1
            for f1 = 1:nf
                fOnset{f1} = job.onset_files{f1,1};
                fOnsetOK = 1;
            end
            doPermute = 1;
        else
            disp('Too many subjects -- incorrect usage -- onsets should be specified earlier');
            out = [];
        end
    end
end
if doPermute
    for Idx=1:size(job.NIRSmat,1)
        %Load NIRS.mat information
        try
            NIRS = [];
            load(job.NIRSmat{Idx,1});
            if ~fOnsetOK
                for f1=1:length(NIRS.Dt.fir.rons)
                    fOnset{f1} = NIRS.Dt.fir.rons{f1};
                end
            else
                
            end
            permute_onsets = 0;
            load(fOnset{1});
            names2 = names;
            GroupOnsets = 0;
            for j=2:nf
                load(fOnset{j});
                if size(names2,2) == size(names,2)
                    for m=1:size(names,2)
                        if ~strcmpi(names2{m},names{m})
                            permute_onsets = 1;
                        end
                    end
                else
                    disp(strvcat(['Different number of onsets in session ' int2str(j)],...
                        'Onsets will be grouped into the same type'));
                    GroupOnsets = 1;
                    %Don't bother trying to permute onsets if number of onsets
                    %doesn't match
                    permute_onsets = 0;
                end
            end
            
            if permute_onsets && ~GroupOnsets
                %Onset file, Movement
                load(fOnset{1});
                names2 = names;
                clear onsets durations
                %generate all permutations
                perm1 = perms(1:size(names2,2));
                %loop over sessions
                for j=2:nf
                    load(fOnset{j});
                    %loop over permutations
                    for k=1:size(perm1,1)
                        %start assuming this is the good permutation
                        good_perm = 1;
                        %loop over onset types
                        for m=1:size(names2,2)
                            %find the permutation for which all the onset names
                            %match the ones of the first session
                            if ~strcmpi(names{perm1(k,m)},names2{m})
                                good_perm = 0; %not a good permutation
                            end
                        end
                        if good_perm == 1, break; end
                    end
                    %If we get here, that's because for this value of k,
                    %there was no mismatch of names, so this is the correct
                    %permutation
                    names = names2;
                    for m=1:size(names2,2)
                        temp_onsets{m} = onsets{perm1(k,m)};
                        temp_durations{m} = durations{perm1(k,m)};
                    end
                    onsets = temp_onsets;
                    durations = temp_durations;

                end
            end
        catch
            disp('Problem permuting files');
        end
    end
end
out.NIRSmat = job.NIRSmat;
end