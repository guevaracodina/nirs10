function out = nirs_run_runVOIRE(job)
%__________________________________________________________________
% Copyright (C) 2010 Laboratoire d'Imagerie Optique et Moleculaire
% Clément Bonnéry
% 2010-12

% le rythme cardiaque me donne
%Friston_nonlinear_2000 :
% rCBF lineaire avec synaptic activity
% CBV grossièrement égal à HbT : différence est dans la plasme qui est
% l'autre partie du sang !!
% Cox : oxygénation du sang : grossièrement HbO/HbT (on tient pas compte du plasma...)
% Débit cardiaque = freq card * volume éjecté
prefix = 'h';

if job.heart_pace==1
    outHP = nirs_run_criugm_paces(job.criugm_paces1);
end

load(job.NIRSmat{:});

jobHb.age = NIRS.Dt.s.age;
jobHb.NIRSmat = job.NIRSmat;
jobHb.Normalize_OD = 0;
jobHb.subject_age = NIRS.Dt.s.age;
jobHb.PVF = [50;50];
jobHb.threshold = 0.1;
nirs_lpf2.lpf_gauss2.fwhm1 =1.5;
nirs_lpf2.lpf_gauss2.downsamplingFactor =1;
nirs_lpf2.lpf_gauss2.downsampleWhen = 1;
jobHb.nirs_lpf2 = nirs_lpf2;

outHb = nirs_run_ODtoHbOHbR(jobHb);

clear NIRS;
load(job.NIRSmat{:});
%use last step of preprocessing
lst = length(NIRS.Dt.fir.pp);
rDtp = NIRS.Dt.fir.pp(lst).p; % path for files to be processed
try
    %loop over data files
    for f=1:size(rDtp,1)
        c = fopen_NIR(rDtp{f});
        [dummy,fil1,dummy2] = fileparts(rDtp{f});
        
        % Attention deux pb
        %-- on a parfois une mauvaise detection des rythmes cardiaques et
        %ceux-ci sont consideres comme bons...
        %-- il faut raisonner en temporel pour voir quels rythmes se
        %correspondent d'une paire a l'autre et ainsi moyenner que ce qui
        %est bon
        
        %>> on cherche une paire qui est fiable (critere : elle n'a pas de grands sauts...)
        %--- on cherche le fichier du rythme qui correspond au fichier de
        %donnees
        try
            fi=1;
            while isempty(strfind(NIRS.Dt.fir.ht{fi,1}.from,fil1(3:end)))
                fi =fi+1;
            end
            hp = NIRS.Dt.fir.ht{fi,1}.pace;
            %--- on charge le rythme cardiaque et on cherche une paire fiable
            %(une evolution sans coupure)
            nc = size(hp,1);
            ci=1;
            cref=0;
            err =[];
            while ci<nc && cref==0
                if sum(hp(ci,:)~= zeros(1,size(hp,2)))>100 && cref==0
                    for t=1:size(hp,2)-1
                        if hp(ci,t+1)>hp(ci,t)+5 || hp(ci,t+1)<hp(ci,t)-5% && t<size(hp,2)-1;
                            err = [err,1];
                        end
                    end
                    if size(err,2)<10
                        cref=ci;
                    end
                    
                end
                ci = ci+1;
                err =[];
            end
            
            if cref==0
                disp([rDtp{f} ' : no heart pace in any channel !'])
            else
                %--- on compare les autres paires et on ne somme que les moments
                %des paires qui ont un rythme cardiaque qui correspond a la paire
                %de reference
                
                %tableau de bouleens qui contient tous les instants a sommer
                cok_t = zeros(size(hp));
                for ci =1:nc
                    for t=1:size(hp,2)
                        if hp(ci,t)< hp(cref,t)+5 && hp(ci,t)> hp(cref,t)-5
                            cok_t(ci,t) =1;
                        end
                    end
                end
                
                Cid = NIRS.Cf.H.C.id;
                COx = zeros(size(c,1)/2,size(c,2));
                %             test avec _cok : voir si en n'ajoutant que les parties de
                %             signal qui sont bonnes
                %             COx_Cok = zeros(size(c,1)/2,size(c,2));
                %             Cok_t = [cok_t;cok_t];
                %             c_Cok = c.*Cok_t;
                for iC = 1:size(Cid,2)/2
                    COx(iC,:) = c(iC,:)./(c(iC,:)+c(iC+size(c,1)/2,:));
                    %                 COx_Cok(iC,:) = c_Cok(iC,:)./(c_Cok(iC,:)+c_Cok(iC+size(c,1)/2,:));
                end
                COx(isnan(COx))=0;
                %             COx_Cok(isnan(COx_Cok))=0;
                %             for i=1:size(COx_Cok,1)
                %                 COx_wCok(i,:) = COx_Cok(i,:)/norm(COx_Cok(i,:));
                %             end
                %             COx_wCok(isnan(COx_wCok))=0;
                
                %test 2 : on ne conserve que les paires valables tout au long
                test = sum(cok_t,2);
                count=0;
                COxsum = zeros(1,size(hp,2));
                for i=1:size(test,1)/2
                    %             if test(i,1)>size(hp,2)-1000
                    COxsum = COxsum +COx(i,:);
                    count = count+1;
                    %             end
                end
                COxsum = COxsum/count;
                
%                 save(fullfile(NIRS.Dt.s.p,['Cox' fil1 '.mat']),'COx');
            end
        catch exception
            disp(exception.identifier);
            disp(exception.stack(1));
            error(msg);
            disp([rDtp{f} ' : heart regressor has not been calculated']);
        end
    end
    
catch
    disp('Could not evaluate COx');
end
out.NIRSmat = job.NIRSmat;