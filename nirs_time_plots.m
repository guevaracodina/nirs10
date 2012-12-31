function nirs_time_plots(d,fs,NC,f,newNIRSlocation,ftype,conc)
if isnumeric(conc{1})
    conc{1} = num2str(conc{1});
    conc{2} = num2str(conc{2});
end

[dirNIRS fil0] = fileparts(newNIRSlocation);
lp = linspace(0,size(d,2)/fs,size(d,2));
Fdir = fullfile(dirNIRS,'dataPlots');
if ~exist(Fdir,'dir'), mkdir(Fdir); end
fname = fullfile(Fdir,[ftype '_' conc{1} '_Sess' int2str(f)]);
private_plot(d(1:NC/2,:),lp,fname)
fname = fullfile(Fdir,[ftype '_' conc{2} '_Sess' int2str(f)]);
private_plot(d(NC/2+1:NC,:),lp,fname)

dt = 30;
lp2 = linspace(0,dt,round(dt*fs));
fname = fullfile(Fdir,[ftype '_' conc{1} '_Sess' int2str(f) '_zoom']);
private_plot(d(1:NC/2,1:round(dt*fs)),lp2,fname)
fname = fullfile(Fdir,[ftype '_' conc{2} '_Sess' int2str(f) '_zoom']);
private_plot(d(NC/2+1:NC,1:round(dt*fs)),lp2,fname)

function private_plot(d,lp,fname)
h = figure; plot(lp,d');
print(h,'-dpng',[fname '.png'],'-r300');
saveas(h,[fname '.fig']);
close(h);