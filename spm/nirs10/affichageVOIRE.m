path = 'D:\Users\Clément\Projets_CRIUGM\Donnees\VOIRE\03\Stat';
subj = '03';
subju ='03';

figure

%% etude 1 : VOIRE
try
    COx_40 =load([path subj '\Coxhr' subju '40.mat']);
    COx_40.COx(1,:) = sum(COx_40.COx,1)/28;
    
    subplot(4,1,2)
    plot(COx_40.COx(1,:));
    ylim([0.75 0.8])
    % ylim([0.75 0.79])
catch
    disp('erreur');
end
try
    COx_60 =load([path subj '\Coxhr' subju '60.mat']);
    COx_60.COx(1,:) = sum(COx_60.COx,1)/28;
    
    subplot(4,1,3)
    plot(COx_60.COx(1,:));
    ylim([0.75 0.8])
    % ylim([0.75 0.79])
catch
    disp('erreur');
end
try
    COx_85 =load([path subj '\Coxhr' subju '80.mat']);
    COx_85.COx(1,:) = sum(COx_85.COx,1)/28;
    
    
    subplot(4,1,4)
    plot(COx_85.COx(1,:));
    ylim([0.75 0.8])
    % ylim([0.75 0.79])
catch
    disp('erreur');
end


%% etude 2 : VITE
try
    COx_re =load([path subj '\Coxhr' subju 'stroopbaseline.mat']);
    
    t = [1:size(COx_re.COx,2)]*0.04;
    COx_re.COx(1,:) = sum(COx_re.COx,1)/28;
    
    subplot(4,1,1)
    plot(t,COx_re.COx(1,:));
    ylim([0.75 0.8])
    % ylim([0.78 0.82])
    
    line([60 60],[0.7 1]);
    line([120 120],[0.7 1]);
    line([180 180],[0.7 1]);
    line([240 240],[0.7 1]);
    line([300 300],[0.7 1]);
    line([360 360],[0.7 1]);
    line([420 420],[0.7 1]);
    line([480 480],[0.7 1]);
catch
    disp('erreur');
end

try
    COx_40 =load([path subj '\Coxhr' subju 'stroop40.mat']);
    COx_40.COx(1,:) = sum(COx_40.COx,1)/28;
    t = [1:size(COx_40.COx,2)]*0.04;
    subplot(4,1,2)
    plot(t,COx_40.COx(1,:));
    ylim([0.75 0.8])
    
    line([60 60],[0.7 1]);
    line([120 120],[0.7 1]);
    line([180 180],[0.7 1]);
    line([240 240],[0.7 1]);
    line([300 300],[0.7 1]);
    line([360 360],[0.7 1]);
    line([420 420],[0.7 1]);
    line([480 480],[0.7 1]);
    % ylim([0.75 0.79])
catch
    disp('erreur');
end
try
    COx_60 =load([path subj '\Coxhr' subju 'stroop60.mat']);
    COx_60.COx(1,:) = sum(COx_60.COx,1)/28;
    
    t = [1:size(COx_60.COx,2)]*0.04;
    subplot(4,1,3)
    plot(t,COx_60.COx(1,:));
    ylim([0.75 0.8])
    line([60 60;120 120;180 180;240 240;300 300;360 360;420 420;480 480],[0.7 1;0.7 1;0.7 1;0.7 1;0.7 1;0.7 1;0.7 1;0.7 1])
    % ylim([0.75 0.79])
    
    line([60 60],[0.7 1]);
    line([120 120],[0.7 1]);
    line([180 180],[0.7 1]);
    line([240 240],[0.7 1]);
    line([300 300],[0.7 1]);
    line([360 360],[0.7 1]);
    line([420 420],[0.7 1]);
    line([480 480],[0.7 1]);
catch
    disp('erreur');
end
try
    COx_85 =load([path subj '\Coxhr' subju 'stroop85.mat']);
    COx_85.COx(1,:) = sum(COx_85.COx,1)/28;
    
    t = [1:size(COx_85.COx,2)]*0.04;
    subplot(4,1,4)
    plot(t,COx_85.COx(1,:));
    ylim([0.75 0.8])
    line([60 60;120 120;180 180;240 240;300 300;360 360;420 420;480 480],[0.7 1;0.7 1;0.7 1;0.7 1;0.7 1;0.7 1;0.7 1;0.7 1])
    % ylim([0.75 0.79])
    
    line([60 60],[0.7 1]);
    line([120 120],[0.7 1]);
    line([180 180],[0.7 1]);
    line([240 240],[0.7 1]);
    line([300 300],[0.7 1]);
    line([360 360],[0.7 1]);
    line([420 420],[0.7 1]);
    line([480 480],[0.7 1]);
catch
    disp('erreur');
end

%% GLM
figure;
subplot(4,1,1);stem(SPM.xXn{1,1}.t(1,:));title(['Sujet ' subj ' - stroop - 1,1' 'Baseline']);
subplot(4,1,2);stem(SPM.xXn{1,1}.t(2,:));title('Denomination');
subplot(4,1,3);stem(SPM.xXn{1,1}.t(3,:));title('Inhib_switch');
subplot(4,1,4);stem(SPM.xXn{1,1}.t(4,:));title('Constante');
figure;
subplot(4,1,1);stem(SPM.xXn{1,2}.t(1,:));title(['Sujet ' subj ' - stroop - 1,2' 'Baseline']);
subplot(4,1,2);stem(SPM.xXn{1,2}.t(2,:));title('Denomination');
subplot(4,1,3);stem(SPM.xXn{1,2}.t(3,:));title('Inhib_switch');
subplot(4,1,4);stem(SPM.xXn{1,2}.t(4,:));title('Constante');
figure;
subplot(4,1,1);stem(SPM.xXn{1,3}.t(1,:));title(['Sujet ' subj ' - stroop - 1,3' 'Baseline']);
subplot(4,1,2);stem(SPM.xXn{1,3}.t(2,:));title('Denomination');
subplot(4,1,3);stem(SPM.xXn{1,3}.t(3,:));title('Inhib_switch');
subplot(4,1,4);stem(SPM.xXn{1,3}.t(4,:));title('Constante');
figure;
subplot(4,1,1);stem(SPM.xXn{1,4}.t(1,:));title(['Sujet ' subj ' - stroop - 1,4' 'Baseline']);
subplot(4,1,2);stem(SPM.xXn{1,4}.t(2,:));title('Denomination');
subplot(4,1,3);stem(SPM.xXn{1,4}.t(3,:));title('Inhib_switch');
subplot(4,1,4);stem(SPM.xXn{1,4}.t(4,:));title('Constante');