function Figure1(A)
%   Figure1(A) generates Figure 1 of the comment to "DNA damage is a 
%   pervasive cause of sequencing errors, directly confounding variant 
%   identification.", Science 355, 752?756 (2017)
%
%   loads tables from git@github.com:broadinstitute/damage.git
%
%   Copyright 2018 Chip Stewart,  The Broad Institute 

X=load_tsv('preAdapterMetrics.noContext.GIV.24Apr2017.tsv')
X=trimStruct(X,strfindk(X.sample_id,'PAAD','v'))
X.QSCORE=-10*log10(X.ERROR_RATE)
k=find(ismember(X.REFALT,'GT'))
k1=find(ismember(X.REFALT,'CA'))

% transform picard to GIV score
top=X.PRO_ALT_BASES./(X.PRO_ALT_BASES+X.PRO_REF_BASES);
bot=X.CON_ALT_BASES./(X.CON_ALT_BASES+X.CON_REF_BASES);

X.picard_GIV=top./bot;

median_null_error_rate=median(bot(k))/1.25
mean_null_error_rate=mean(bot(k))

clf
subplot(2,3,1)
semilogy(X.QSCORE(k),X.GIV0(k),'.')
xlim([20 60])
SETS=sort(unique(X.SET))
SETS=flipud(SETS)
hold on
for i=2:length(SETS)
    k1=find(ismember(X.SET(k),SETS(i)));
    semilogy(X.QSCORE(k(k1)),X.GIV0(k(k1)),'.')
end
ylim([0.8 4e1])
p1=patch([20 30 30 20 20],[2 2 40 40 2],'black','FaceAlpha',0.13,'EdgeColor','none')
co=get(0,'DefaultAxesColorOrder')
for i=2:length(SETS)
    k1=find(ismember(X.SET(k),SETS(i)));
    semilogy(X.QSCORE(k(k1)),X.GIV0(k(k1)),'.','color',co(i,:))
end
hold off
xlabel('oxoQ')
ylabel({'GIV_{G\_T} reported'},'interpreter','tex')
grid on;
[h,icons,plots,legend_text]=legend(SETS)

NSETS=length(SETS)
kl=findobj(icons,'type','Line')
for o=kl(2:2:(2*NSETS))'
    o.MarkerSize=15
end
kt=findobj(icons,'type','Text')
for o=kt'
    o.FontSize=12
end

x=20:60;
y=1+10.^(-x/10)/median_null_error_rate;
text(-0.2,1.05,'A','units','normalized','fontsize',18)


subplot(2,3,2)
semilogy(X.QSCORE(k),X.GIV(k),'.')
xlim([20 60])

hold on
for i=2:length(SETS)
    k1=find(ismember(X.SET(k),SETS(i)));
    semilogy(X.QSCORE(k(k1)),X.GIV(k(k1)),'.')
end
ylim([0.8 4e1])
p1=patch([20 30 30 20 20],[2 2 40 40 2],'black','FaceAlpha',0.13,'EdgeColor','none')
for i=2:length(SETS)
    k1=find(ismember(X.SET(k),SETS(i)));
    semilogy(X.QSCORE(k(k1)),X.GIV(k(k1)),'.','color',co(i,:))
end

hold off
xlabel('oxoQ')
ylabel({'GIV_{G\_T} corrected'},'interpreter','tex')


[h,icons,plots,legend_text]=legend(SETS)

NSETS=length(SETS)
kl=findobj(icons,'type','Line')
for o=kl(2:2:(2*NSETS))'
    o.MarkerSize=15
end
kt=findobj(icons,'type','Text')
for o=kt'
    o.FontSize=12
end

grid on;

text(-0.2,1.05,'B','units','normalized','fontsize',18)


set(gcf,'Position',[10 300 1200 550])

%
subplot(2,3,3)
Q=load('bqHisto.11Apr2017.mat')

[u,k]=sort(Q.oxoQ)
Q=trimStruct(Q,k)

k1=2;
k2=find( abs(Q.oxoQ-48.882)<0.01  & abs(Q.meanbq-33.1532)<0.01)
b1=bar(Q.x(k1).bq,Q.x(k1).count/1e9);
b1.FaceColor = [0 0 1];
b1.FaceAlpha = 0.05;

p1=patch([15 30 30 15 15],[0 0 1.5 1.5 0],'black','FaceAlpha',0.13,'EdgeColor','none')

hold on;
b1=bar(Q.x(k1).bq,Q.x(k1).count/1e9);
b1.FaceColor = [0 0 1];
b1.FaceAlpha = 0.6;



b2=bar(Q.x(k2).bq,Q.x(k2).count/1e9);
b2.FaceColor = [1 0 0];
b2.FaceAlpha = 0.6;
xlim([15 50])
xlabel('base quality')
ylabel('bases (billions)')
text(0.015,0.95,Q.sample_id{k1},'units','normalized','color','b','fontsize',12)
text(0.015,0.87,sprintf('oxoQ=%.1f',Q.oxoQ(k1)),'units','normalized','color','b','fontsize',12)
text(0.57,0.95,Q.sample_id{k2},'units','normalized','color','r','fontsize',12)
text(0.57,0.87,sprintf('oxoQ=%.1f',Q.oxoQ(k2)),'units','normalized','color','r','fontsize',12)
grid on
hold off;
 
text(-0.2,1.05,'C','units','normalized','fontsize',18)

%
subplot(2,3,6)
k=find(ismember(X.REFALT,'GT'))
x=X.GIV0_TOT(k)./X.GIV_TOT(k); x(x>1)=1;
semilogx(X.GIV(k),x,'ko','markersize',3); 
xlabel({'GIV_{G\_T} corrected '},'interpreter','tex')
ylabel('fraction of bases selected')
xlim([0.9 40])
grid on

p1=patch([5 40 40 5 5],[0 0 0.25 0.25 0],'black','FaceAlpha',0.13,'EdgeColor','none')
hold on;
semilogx(X.GIV(k),x,'ko','markersize',3); 
hold off
text(-0.2,1.05,'F','units','normalized','fontsize',18)

i=1
%% make lego input mat file 
for i=1:2
    if (i==1), v='revised'; panel='E'; label ='corrected'; ipanel=5; end
    if (i==2), v='default'; panel='D'; label ='reported'; ipanel=4; end
    if exist(['data_' v '/estimate_damage_location_context_for_R.oxoQ.le.30.8Jan2018.mat'],'file'), continue; end
    D=load(['data_' v '/estimate_damage_location_context_for_R.14Apr2017.mat'])

    
    D.DERROR_RATE =  (D.PRO_ALT_BASES-D.CON_ALT_BASES)./(D.PRO_ALT_BASES+D.CON_ALT_BASES + D.PRO_REF_BASES+D.CON_REF_BASES);
    D.PRO_ERROR_RATE =  D.PRO_ALT_BASES./(D.PRO_ALT_BASES + D.PRO_REF_BASES);
    D.CON_ERROR_RATE =  D.CON_ALT_BASES./(D.CON_ALT_BASES + D.CON_REF_BASES);
    D.GIV =  D.PRO_ERROR_RATE./D.CON_ERROR_RATE;
    
    de=D.DERROR_RATE;
    D.QSCORE = -10*log10(abs(de));
    k=find(de<0);
    D.QSCORE(k) = 100; %-D.QSCORE(k);
    
    k=find(ismember(D.sample,X.sample_id));
    D=trimStruct(D,k)

    % over all 1235
    S0=unique(X.sample_id)
    D=trimStruct(D,ismember(D.sample,S0))
    
    % OXOG <30
    k=find(ismember(X.REF_BASE,'G').*ismember(X.ALT_BASE,'T').*(abs(X.QSCORE)<=30))
    X1=trimStruct(X,k)
    k=find(ismember(D.sample,X1.sample_id))
    D1=trimStruct(D,k)
    if (0)
        kx=find(ismember(D1.REF_BASE,'G').*ismember(D1.ALT_BASE,'T'))
        DX=trimStruct(D1,kx)
        DX.QSCORE
    end
    save(['data_' v '/estimate_damage_location_context_for_R.oxoQ.le.30.8Jan2018.mat'],'-struct','D1') 
        
    if(i==1),D0=D1; end
end

%%
P=[]
P.PROCON='CON';


v='default'; panel='D'; label ='reported'; ipanel=4;
subplot(2,3,ipanel)

D=load(['data_' v '/estimate_damage_location_context_for_R.oxoQ.le.30.8Jan2018.mat'])
% over all 1235 
[XY,XY1,hb,hl,hp,h1]=plotLegoFromDetailMetrics1(D,P)
text(hb,-0.1,1.05,'D','units','normalized','fontsize',18)

v='revised'; panel='E'; label ='corrected'; ipanel=5;
subplot(2,3,ipanel)

D=load(['data_' v '/estimate_damage_location_context_for_R.oxoQ.le.30.8Jan2018.mat'])
if (0)
    kx=find(ismember(D.REF_BASE,'G').*ismember(D.ALT_BASE,'T'))
    DX1=trimStruct(D,kx)
    DX1.QSCORE
end

% over all 1235 
[XY,XY1,hb,hl,hp,h1]=plotLegoFromDetailMetrics1(D,P)
text(hb,-0.1,1.05,'E','units','normalized','fontsize',18)

fplt=['plots/Fig1.' TODAY ]
saveas(gcf,[fplt '.png'],'png')
print([fplt '.eps'],'-depsc','-painters')
print([fplt '.svg'],'-dsvg','-painters')


%%
    

%%  make table
REFALT=unique(X.REFALT)
k=find(cellfun(@length,regexp(REFALT,'.T|.G','match'))>0)
REFALT=REFALT(k)
X1=trimStruct(X,ismember(X.REFALT,'GT'))

REFALTLAB={'G>T','C>T','T>A','G>C','A>C' 'T>C'}
REFALT=regexprep(REFALTLAB,'>','')
X=trimStruct(X,ismember(X.REFALT,REFALT))
SETS=sort(unique(X.SET))
T=load_tsv('GDC.TCGA.bams.name.sequencing_date.24Apr2017.tsv')
T.id0=T.id;
T.proj=mysubstring(T.id,1,4);
T.id=mysubstring(T.id,6,10);
T.id=regexprep(T.id,'-01$','-TP')
T.id=regexprep(T.id,'-11$','-NT')
T.id=regexprep(T.id,'-10$','-NB')
T.id=regexprep(T.id,'-06$','-TM')
T.id=regexprep(T.id,'-03$','-TB')
T.short_sample_code=mysubstring(T.id,9,2)
T=trimStruct(T,ismember(T.short_sample_code,{'TB','TM','TP'}))
[i m]=ismember(X.id,T.id)
%tab(i)
X.GDC_T1=NaN*X.GIV;
X.GDC_T2=NaN*X.GIV;
X.GDC_T1(i)=T.t1(m(i));
X.GDC_T2(i)=T.t2(m(i));
plot(X.bam_start_date,X.GDC_T1,'+')
X1=trimStruct(X,ismember(X.REFALT,'GT'))
X1=rmfield(X1,{'CONTEXT','LIBRARY'})
X1.BAM_START_DATE=repmat({''},size(X1.SAMPLE_ALIAS))
k=find(X1.GDC_T1>0)
X1.BAM_START_DATE(k)=cellstr(datestr(X1.GDC_T1(k), 'ddmmmyyyy'))
X1.BAM_END_DATE=repmat({''},size(X1.SAMPLE_ALIAS))
k=find(X1.GDC_T2>0)
X1.BAM_END_DATE(k)=cellstr(datestr(X1.GDC_T2(k), 'ddmmmyyyy'))

for i=1:length(REFALT)
    REFALT1=REFALT{i}
    x1=trimStruct(X,ismember(X.REFALT,REFALT1))
    F1=['PRO_REF_BASES_' REFALT1]    
    X1.(F1)=x1.PRO_REF_BASES;
    F1=['PRO_ALT_BASES_' REFALT1]    
    X1.(F1)=x1.PRO_ALT_BASES;
    F1=['CON_REF_BASES_' REFALT1]    
    X1.(F1)=x1.CON_REF_BASES;
    F1=['CON_ALT_BASES_' REFALT1]    
    X1.(F1)=x1.CON_ALT_BASES;
    F1=['ERROR_RATE_' REFALT1]    
    X1.(F1)=x1.ERROR_RATE;
    F1=['ARTQSCORE_' REFALT1]    
    X1.(F1)=x1.QSCORE;
    F1=['PICARD_GIV_' REFALT1]
    X1.(F1)=x1.picard_GIV;
    F1=['GIV_' REFALT1 '_DEFAULT']
    X1.(F1)=x1.GIV0;
    F1=['GIV_' REFALT1 '_QUALITYSCORE20']
    X1.(F1)=x1.GIV;
end
X1=rmfield(X1,{'REF_BASE','ALT_BASE','PRO_REF_BASES','PRO_ALT_BASES','CON_REF_BASES','CON_ALT_BASES','ERROR_RATE','QSCORE','sample_id','id','REFALT','GIV','GIV0','GIV0_TOT','GIV_TOT','GIV_FAM','picard_GIV'})
X1=rmfield(X1,{'bam_start_date','bam_end_date','GDC_T1','GDC_T2'})
printStruct(X1,-1,['TableS1.Picard.GIV.' TODAY '.tsv'])
