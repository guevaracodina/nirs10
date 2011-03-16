%% Affichage rythme cardiaque

figure;
subplot(4,1,1);
t = ((1:size(NIRS.Dt.fir.Sess(1,1).fR{1,1},1))*0.04)';
plot(t,NIRS.Dt.fir.Sess(1,1).fR{1,1});
subplot(4,1,2);
t = ((1:size(NIRS.Dt.fir.Sess(1,2).fR{1,1},1))*0.04)';
plot(t,NIRS.Dt.fir.Sess(1,2).fR{1,1});
subplot(4,1,3);
t = ((1:size(NIRS.Dt.fir.Sess(1,3).fR{1,1},1))*0.04)';
plot(t,NIRS.Dt.fir.Sess(1,3).fR{1,1});
subplot(4,1,4);
t = ((1:size(NIRS.Dt.fir.Sess(1,4).fR{1,1},1))*0.04)';
plot(t,NIRS.Dt.fir.Sess(1,4).fR{1,1});
title('nath');

%% Affichage des moyennes de HbO sur le temps
p ='D:\Users\Clément\Projets_CRIUGM\Donnees\Said_etude1\aurelia';

load(fullfile(p,'NIRS.mat'));
figure;
load(fullfile(p,'hbraurelia.nirs'),'-mat');
Cok = NIRS.Cf.H.C.ok;
t = (1:size(d,1))*0.04;
som=zeros(size(d,1),1);
for i=1:size(Cok,2)/2
    som = som+d(:,Cok(1,i));
end
som = som/(size(Cok,2)/2);
plot(t,som);title('AURELIA');

%%
p ='D:\Users\Clément\Projets_CRIUGM\Donnees\Said_etude1\catherine';

load(fullfile(p,'NIRS.mat'));
figure;
load(fullfile(p,'hbrcath.nirs'),'-mat');
Cok = NIRS.Cf.H.C.ok;
t = (1:size(d,1))*0.04;
som=zeros(size(d,1),1);
for i=1:size(Cok,2)/2
    som = som+d(:,Cok(1,i));
end
som = som/(size(Cok,2)/2);
subplot(4,1,1);plot(t,som);hold on; 

load(fullfile(p,'hbrcath40.nirs'),'-mat');
Cok = NIRS.Cf.H.C.ok;
t = (1:size(d,1))*0.04;
som=zeros(size(d,1),1);
for i=1:size(Cok,2)/2
    som = som+d(:,Cok(1,i));
end
som = som/(size(Cok,2)/2);
subplot(4,1,2);plot(t,som);hold on;

load(fullfile(p,'hbrcath60.nirs'),'-mat');
Cok = NIRS.Cf.H.C.ok;
t = (1:size(d,1))*0.04;
som=zeros(size(d,1),1);
for i=1:size(Cok,2)/2
    som = som+d(:,Cok(1,i));
end
som = som/(size(Cok,2)/2);
subplot(4,1,3);plot(t,som);hold on;

load(fullfile(p,'hbrcath85.nirs'),'-mat');
Cok = NIRS.Cf.H.C.ok;
t = (1:size(d,1))*0.04;
som=zeros(size(d,1),1);
for i=1:size(Cok,2)/2
    som = som+d(:,Cok(1,i));
end
som = som/(size(Cok,2)/2);
subplot(4,1,4);plot(t,som);title('CATHERINE');

%%
p ='D:\Users\Clément\Projets_CRIUGM\Donnees\Said_etude1\claudine';

load(fullfile(p,'NIRS.mat'));
figure;
load(fullfile(p,'hbrclaudineHmax.nirs'),'-mat');
Cok = NIRS.Cf.H.C.ok;
t = (1:size(d,1))*0.04;
som=zeros(size(d,1),1);
for i=1:size(Cok,2)/2
    som = som+d(:,Cok(1,i));
end
som = som/(size(Cok,2)/2);
subplot(4,1,1);plot(t,som);hold on;

load(fullfile(p,'hbrClaudineH40.nirs'),'-mat');
Cok = NIRS.Cf.H.C.ok;
t = (1:size(d,1))*0.04;
som=zeros(size(d,1),1);
for i=1:size(Cok,2)/2
    som = som+d(:,Cok(1,i));
end
som = som/(size(Cok,2)/2);
subplot(4,1,2);plot(t,som);hold on;

load(fullfile(p,'hbrclaudineH60.nirs'),'-mat');
Cok = NIRS.Cf.H.C.ok;
t = (1:size(d,1))*0.04;
som=zeros(size(d,1),1);
for i=1:size(Cok,2)/2
    som = som+d(:,Cok(1,i));
end
som = som/(size(Cok,2)/2);
subplot(4,1,3);plot(t,som);hold on;

load(fullfile(p,'hbrcLAUDINEh85.nirs'),'-mat');
Cok = NIRS.Cf.H.C.ok;
t = (1:size(d,1))*0.04;
som=zeros(size(d,1),1);
for i=1:size(Cok,2)/2
    som = som+d(:,Cok(1,i));
end
som = som/(size(Cok,2)/2);
subplot(4,1,4);plot(t,som);title('CLAUDINE');

%%
p ='D:\Users\Clément\Projets_CRIUGM\Donnees\Said_etude1\clement';

load(fullfile(p,'NIRS.mat'));
figure;
load(fullfile(p,'hbrclement1.nirs'),'-mat');
Cok = NIRS.Cf.H.C.ok;
t = (1:size(d,1))*0.04;
som=zeros(size(d,1),1);
for i=1:size(Cok,2)/2
    som = som+d(:,Cok(1,i));
end
som = som/(size(Cok,2)/2);
subplot(4,1,1);plot(t,som);hold on;

load(fullfile(p,'hbrCB40.nirs'),'-mat');
Cok = NIRS.Cf.H.C.ok;
t = (1:size(d,1))*0.04;
som=zeros(size(d,1),1);
for i=1:size(Cok,2)/2
    som = som+d(:,Cok(1,i));
end
som = som/(size(Cok,2)/2);
subplot(4,1,2);plot(t,som);hold on;

load(fullfile(p,'hbrCB60.nirs'),'-mat');
Cok = NIRS.Cf.H.C.ok;
t = (1:size(d,1))*0.04;
som=zeros(size(d,1),1);
for i=1:size(Cok,2)/2
    som = som+d(:,Cok(1,i));
end
som = som/(size(Cok,2)/2);
subplot(4,1,3);plot(t,som);hold on;

load(fullfile(p,'hbrCB85.nirs'),'-mat');
Cok = NIRS.Cf.H.C.ok;
t = (1:size(d,1))*0.04;
som=zeros(size(d,1),1);
for i=1:size(Cok,2)/2
    som = som+d(:,Cok(1,i));
end
som = som/(size(Cok,2)/2);
subplot(4,1,4);plot(t,som);title('CLEMENT');
%%
p ='D:\Users\Clément\Projets_CRIUGM\Donnees\Said_etude1\clement';

load(fullfile(p,'NIRS.mat'));
figure;
load(fullfile(p,'hbrCB60.nirs'),'-mat');
Cok = NIRS.Cf.H.C.ok;
t = (1:size(d,1))*0.04;
som=zeros(size(d,1),1);
for i=1:size(Cok,2)/2
    som = som+d(:,Cok(1,i));
end
som = som/(size(Cok,2)/2);
plot(t,som,'r');hold on;

som=zeros(size(d,1),1);
for i=1:size(Cok,2)/2
    som = som+d(:,Cok(1,i)-size(Cok,2)/2);
end
som = som/(size(Cok,2)/2);
plot(t,som);


%%
p ='D:\Users\Clément\Projets_CRIUGM\Donnees\Said_etude1\gb';

load(fullfile(p,'NIRS.mat'));
figure;
load(fullfile(p,'hbrGBMAX1.nirs'),'-mat');
Cok = NIRS.Cf.H.C.ok;
t = (1:size(d,1))*0.04;
som=zeros(size(d,1),1);
for i=1:size(Cok,2)/2
    som = som+d(:,Cok(1,i));
end
som = som/(size(Cok,2)/2);
subplot(4,1,1);plot(t,som);hold on;

load(fullfile(p,'hbrGB40.nirs'),'-mat');
Cok = NIRS.Cf.H.C.ok;
t = (1:size(d,1))*0.04;
som=zeros(size(d,1),1);
for i=1:size(Cok,2)/2
    som = som+d(:,Cok(1,i));
end
som = som/(size(Cok,2)/2);
subplot(4,1,2);plot(t,som);hold on;

load(fullfile(p,'hbrGB60.nirs'),'-mat');
Cok = NIRS.Cf.H.C.ok;
t = (1:size(d,1))*0.04;
som=zeros(size(d,1),1);
for i=1:size(Cok,2)/2
    som = som+d(:,Cok(1,i));
end
som = som/(size(Cok,2)/2);
subplot(4,1,3);plot(t,som);hold on;

load(fullfile(p,'hbrGB85.nirs'),'-mat');
Cok = NIRS.Cf.H.C.ok;
t = (1:size(d,1))*0.04;
som=zeros(size(d,1),1);
for i=1:size(Cok,2)/2
    som = som+d(:,Cok(1,i));
end
som = som/(size(Cok,2)/2);
subplot(4,1,4);plot(t,som);title('GB');

%%
p ='D:\Users\Clément\Projets_CRIUGM\Donnees\Said_etude1\louiseric';

load(fullfile(p,'NIRS.mat'));
figure;
load(fullfile(p,'hbrLouis-Eric.nirs'),'-mat');
Cok = NIRS.Cf.H.C.ok;
t = (1:size(d,1))*0.04;
som=zeros(size(d,1),1);
for i=1:size(Cok,2)/2
    som = som+d(:,Cok(1,i));
end
som = som/(size(Cok,2)/2);
subplot(4,1,1);plot(t,som);hold on;

load(fullfile(p,'hbrLE40.nirs'),'-mat');
Cok = NIRS.Cf.H.C.ok;
t = (1:size(d,1))*0.04;
som=zeros(size(d,1),1);
for i=1:size(Cok,2)/2
    som = som+d(:,Cok(1,i));
end
som = som/(size(Cok,2)/2);
subplot(4,1,2);plot(t,som);hold on;

load(fullfile(p,'hbrLE60.nirs'),'-mat');
Cok = NIRS.Cf.H.C.ok;
t = (1:size(d,1))*0.04;
som=zeros(size(d,1),1);
for i=1:size(Cok,2)/2
    som = som+d(:,Cok(1,i));
end
som = som/(size(Cok,2)/2);
subplot(4,1,3);plot(t,som);hold on;

load(fullfile(p,'hbrLE85.nirs'),'-mat');
Cok = NIRS.Cf.H.C.ok;
t = (1:size(d,1))*0.04;
som=zeros(size(d,1),1);
for i=1:size(Cok,2)/2
    som = som+d(:,Cok(1,i));
end
som = som/(size(Cok,2)/2);
subplot(4,1,4);plot(t,som);title('LOUIS-ERIC');

%%
p ='D:\Users\Clément\Projets_CRIUGM\Donnees\Said_etude1\marjo';

load(fullfile(p,'NIRS.mat'));
figure;
load(fullfile(p,'hbrmarjo.nirs'),'-mat');
Cok = NIRS.Cf.H.C.ok;
t = (1:size(d,1))*0.04;
som=zeros(size(d,1),1);
for i=1:size(Cok,2)/2
    som = som+d(:,Cok(1,i));
end
som = som/(size(Cok,2)/2);
subplot(4,1,1);plot(t,som);hold on;

load(fullfile(p,'hbrmarjo40.nirs'),'-mat');
Cok = NIRS.Cf.H.C.ok;
t = (1:size(d,1))*0.04;
som=zeros(size(d,1),1);
for i=1:size(Cok,2)/2
    som = som+d(:,Cok(1,i));
end
som = som/(size(Cok,2)/2);
subplot(4,1,2);plot(t,som);hold on;

load(fullfile(p,'hbrmarjo60.nirs'),'-mat');
Cok = NIRS.Cf.H.C.ok;
t = (1:size(d,1))*0.04;
som=zeros(size(d,1),1);
for i=1:size(Cok,2)/2
    som = som+d(:,Cok(1,i));
end
som = som/(size(Cok,2)/2);
subplot(4,1,3);plot(t,som);hold on;

load(fullfile(p,'hbrmarjo85.nirs'),'-mat');
Cok = NIRS.Cf.H.C.ok;
t = (1:size(d,1))*0.04;
som=zeros(size(d,1),1);
for i=1:size(Cok,2)/2
    som = som+d(:,Cok(1,i));
end
som = som/(size(Cok,2)/2);
subplot(4,1,4);plot(t,som);title('MARJO');

%%
p ='D:\Users\Clément\Projets_CRIUGM\Donnees\Said_etude1\maxime';

load(fullfile(p,'NIRS.mat'));
figure;
load(fullfile(p,'hbrmax1.nirs'),'-mat');
Cok = NIRS.Cf.H.C.ok;
t = (1:size(d,1))*0.04;
som=zeros(size(d,1),1);
for i=1:size(Cok,2)/2
    som = som+d(:,Cok(1,i));
end
som = som/(size(Cok,2)/2);
subplot(4,1,1);plot(t,som);hold on;

load(fullfile(p,'hbrMax40.nirs'),'-mat');
Cok = NIRS.Cf.H.C.ok;
t = (1:size(d,1))*0.04;
som=zeros(size(d,1),1);
for i=1:size(Cok,2)/2
    som = som+d(:,Cok(1,i));
end
som = som/(size(Cok,2)/2);
subplot(4,1,2);plot(t,som);hold on;

load(fullfile(p,'hbrMax60.nirs'),'-mat');
Cok = NIRS.Cf.H.C.ok;
t = (1:size(d,1))*0.04;
som=zeros(size(d,1),1);
for i=1:size(Cok,2)/2
    som = som+d(:,Cok(1,i));
end
som = som/(size(Cok,2)/2);
subplot(4,1,3);plot(t,som);hold on;

load(fullfile(p,'hbrMax80.nirs'),'-mat');
Cok = NIRS.Cf.H.C.ok;
t = (1:size(d,1))*0.04;
som=zeros(size(d,1),1);
for i=1:size(Cok,2)/2
    som = som+d(:,Cok(1,i));
end
som = som/(size(Cok,2)/2);
subplot(4,1,4);plot(t,som);title('MAXIME');

%%
p ='D:\Users\Clément\Projets_CRIUGM\Donnees\Said_etude1\michele';

load(fullfile(p,'NIRS.mat'));
figure;
load(fullfile(p,'hbrMichelemax1.nirs'),'-mat');
Cok = NIRS.Cf.H.C.ok;
t = (1:size(d,1))*0.04;
som=zeros(size(d,1),1);
for i=1:size(Cok,2)/2
    som = som+d(:,Cok(1,i));
end
som = som/(size(Cok,2)/2);
subplot(4,1,1);plot(t,som);hold on;

load(fullfile(p,'hbrMichele40.nirs'),'-mat');
Cok = NIRS.Cf.H.C.ok;
t = (1:size(d,1))*0.04;
som=zeros(size(d,1),1);
for i=1:size(Cok,2)/2
    som = som+d(:,Cok(1,i));
end
som = som/(size(Cok,2)/2);
subplot(4,1,2);plot(t,som);hold on;

load(fullfile(p,'hbrMichele60.nirs'),'-mat');
Cok = NIRS.Cf.H.C.ok;
t = (1:size(d,1))*0.04;
som=zeros(size(d,1),1);
for i=1:size(Cok,2)/2
    som = som+d(:,Cok(1,i));
end
som = som/(size(Cok,2)/2);
subplot(4,1,3);plot(t,som);hold on;

load(fullfile(p,'hbrMichele85.nirs'),'-mat');
Cok = NIRS.Cf.H.C.ok;
t = (1:size(d,1))*0.04;
som=zeros(size(d,1),1);
for i=1:size(Cok,2)/2
    som = som+d(:,Cok(1,i));
end
som = som/(size(Cok,2)/2);
subplot(4,1,4);plot(t,som);title('MICHELE');
%%
p ='D:\Users\Clément\Projets_CRIUGM\Donnees\Said_etude1\michele';

load(fullfile(p,'NIRS.mat'));
figure;
load(fullfile(p,'hbrMichele60.nirs'),'-mat');
Cok = NIRS.Cf.H.C.ok;
t = (1:size(d,1))*0.04;
som=zeros(size(d,1),1);
for i=1:size(Cok,2)/2
    som = som+d(:,Cok(1,i));
end
som = som/(size(Cok,2)/2);
plot(t,som,'r');hold on;

som=zeros(size(d,1),1);
for i=1:size(Cok,2)/2
    som = som+d(:,Cok(1,i)-size(Cok,2)/2);
end
som = som/(size(Cok,2)/2);
plot(t,som);


%%
p ='D:\Users\Clément\Projets_CRIUGM\Donnees\Said_etude1\nathalie';

load(fullfile(p,'NIRS.mat'));
figure;
load(fullfile(p,'hbrnathalie.nirs'),'-mat');
Cok = NIRS.Cf.H.C.ok;
t = (1:size(d,1))*0.04;
som=zeros(size(d,1),1);
for i=1:size(Cok,2)/2
    som = som+d(:,Cok(1,i));
end
som = som/(size(Cok,2)/2);
subplot(4,1,1);plot(t,som);hold on;

load(fullfile(p,'hbrNath40.nirs'),'-mat');
Cok = NIRS.Cf.H.C.ok;
t = (1:size(d,1))*0.04;
som=zeros(size(d,1),1);
for i=1:size(Cok,2)/2
    som = som+d(:,Cok(1,i));
end
som = som/(size(Cok,2)/2);
subplot(4,1,2);plot(t,som);hold on;

load(fullfile(p,'hbrNath60.nirs'),'-mat');
Cok = NIRS.Cf.H.C.ok;
t = (1:size(d,1))*0.04;
som=zeros(size(d,1),1);
for i=1:size(Cok,2)/2
    som = som+d(:,Cok(1,i));
end
som = som/(size(Cok,2)/2);
subplot(4,1,3);plot(t,som);hold on;

load(fullfile(p,'hbrNath80.nirs'),'-mat');
Cok = NIRS.Cf.H.C.ok;
t = (1:size(d,1))*0.04;
som=zeros(size(d,1),1);
for i=1:size(Cok,2)/2
    som = som+d(:,Cok(1,i));
end
som = som/(size(Cok,2)/2);
subplot(4,1,4);plot(t,som);title('NATHALIE');

%%
p ='D:\Users\Clément\Projets_CRIUGM\Donnees\Said_etude1\paulolivier';

load(fullfile(p,'NIRS.mat'));
figure;
load(fullfile(p,'hbrPO.nirs'),'-mat');
Cok = NIRS.Cf.H.C.ok;
t = (1:size(d,1))*0.04;
som=zeros(size(d,1),1);
for i=1:size(Cok,2)/2
    som = som+d(:,Cok(1,i));
end
som = som/(size(Cok,2)/2);
subplot(4,1,1);plot(t,som);hold on;

load(fullfile(p,'hbrpo40.nirs'),'-mat');
Cok = NIRS.Cf.H.C.ok;
t = (1:size(d,1))*0.04;
som=zeros(size(d,1),1);
for i=1:size(Cok,2)/2
    som = som+d(:,Cok(1,i));
end
som = som/(size(Cok,2)/2);
subplot(4,1,2);plot(t,som);hold on;

load(fullfile(p,'hbrPO60.nirs'),'-mat');
Cok = NIRS.Cf.H.C.ok;
t = (1:size(d,1))*0.04;
som=zeros(size(d,1),1);
for i=1:size(Cok,2)/2
    som = som+d(:,Cok(1,i));
end
som = som/(size(Cok,2)/2);
subplot(4,1,3);plot(t,som);hold on;

load(fullfile(p,'hbrPO85.nirs'),'-mat');
Cok = NIRS.Cf.H.C.ok;
t = (1:size(d,1))*0.04;
som=zeros(size(d,1),1);
for i=1:size(Cok,2)/2
    som = som+d(:,Cok(1,i));
end
som = som/(size(Cok,2)/2);
subplot(4,1,4);plot(t,som);title('PAUL-OLIVIER');