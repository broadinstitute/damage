function [XY,X1,hb,hl,hp,h1]=plotLegoFromDetailMetrics1(X,P)
%
% plot and save mutation spectrum lego plot for one sample based on maf file
% and coverasge files in standard firehose location
%
% inputs:   PreAdapter DETAIL metrics file: filename or structure with fields:
%               .Reference_Allele
%               .Tumor_Seq_Allele2
%               .ref_context
%               .i_t_ALT_F1R2
%               .PRO_ALT_BASES   (eg. F1ref_CON for G>T)
%               .CON_REF_BASES
%               .CON_ALT_BASES   (eg. F2ref_PRO for C>A)
%           P: optional parameters
%              .zscale = true (print zscale on main rate plot)
%              .label = label for lego to override the first Tumor_Sample_Barcode
%
%
% outputs:  XY: structure for lego plot
%               XY.n:  8x12 matrix of mutation counts per lego bin
%               XY.ncov: 8x12 matrix of base coverage per lego bin
%               XY.cat:  8x12 cell matrix of context+mutation categories per lego bin
%					     eg. G-C>A-T is a C>A mutation after a G before a T
%               XY.col:  8x12x3 matrix of RGB colors per lego bin'
%
% CS:  25Feb2017, based on plotMutationSpectrumCategLegos

fprintf('plotMutationSpectrumCategLegos\n');

if ~exist('P','var'), P=[]; end
P = impose_default_value(P,'zscale',true);
P = impose_default_value(P,'subplot',false);
P = impose_default_value(P,'log',false);
P = impose_default_value(P,'baseline',0);
P = impose_default_value(P,'ztick',[]);
P = impose_default_value(P,'zmax',[]);
P = impose_default_value(P,'PROCON','BOTH');

if nargin<2
    P=[]
end

XFILE='';
if ~isstruct(X) && exist(X,'file')
    fprintf('PreAdapterMetrics file:\t %s\n',X);
    XFILE=X;
    X=load_tsv(X);
else
    fprintf('PreAdapterMetrics structure length:\t %d\n',length(X.CONTEXT));
end
if (length(X.CONTEXT)<1)
    XY.n=[];
    XY.ncov=[];
    XY.cat=[];
    XY.col=[];
    X=MAF;
    text(0.5,0.5,'Nothing to put on lego plot','units','normalized')
    set(gca,'visible','off')
    return;
end

% C categ struct
base='ACGT';
n=0; C=[];

% lego convention A and C REF
FT=[char(X.REF_BASE) char(X.ALT_BASE)];
k=find( (FT(:,1)=='A')|(FT(:,1)=='C'));
Xa=trimStruct(X,k);
FT=FT(k,:);
k=find(FT(:,1)~=FT(:,2));
Xa=trimStruct(Xa,k);
FT=FT(k,:);
CTX=char(Xa.CONTEXT);
Xa.SNP=strcat(Xa.REF_BASE,'>',Xa.ALT_BASE);
Xa.BRAB=strcat(mysubstring(Xa.CONTEXT,1,1),'-',Xa.REF_BASE,'>',Xa.ALT_BASE,'-',mysubstring(Xa.CONTEXT,3,1));
% base ref alt base
BRAB=unique(Xa.BRAB)
N=length(BRAB)
X1=[]
if isfield(X,'SAMPLE_ALIAS')
    s='SAMPLE_ALIAS'
else
    s='sample'
end
NS=length(unique(X.(s)))
for i=1:N
    k=find(ismember(Xa.BRAB,BRAB(i)));
    x1=trimStruct(Xa,k(1));
    x1.PRO_REF_BASES(1) = sum(Xa.PRO_REF_BASES(k));
    x1.PRO_ALT_BASES(1) = sum(Xa.PRO_ALT_BASES(k));
    x1.CON_REF_BASES(1) = sum(Xa.CON_REF_BASES(k));
    x1.CON_ALT_BASES(1) = sum(Xa.CON_ALT_BASES(k));
    x1=trimStruct(x1,1);
    x1.SAMPLE_ALIAS={''};
    x1.LIBRARY={''};
    if isempty(X1)
        X1=x1;
    else
        X1=mergeStruct(X1,x1);
    end
end

X1.ERROR_RATE =  (X1.PRO_ALT_BASES-X1.CON_ALT_BASES)./(X1.PRO_ALT_BASES+x1.CON_ALT_BASES + X1.PRO_REF_BASES+x1.CON_REF_BASES);
e=X1.ERROR_RATE;
X1.QSCORE = -10*log10(abs(e))
k=find(e<0)
X1.QSCORE(k) = 100; %-X1.QSCORE(k)

colors = [1 1 0;0 0.7 0.7;1 0 0;0.1 0.8 0.1;0 0.2 0.8;0.5 0.3 0.7;];

XY=[];
base={'T';'C';'A';'G'};
snp={'C>T';'C>A';'C>G';'A>G';'A>C';'A>T'};
%snp={'G>A';'G>T';'G>C';'T>C';'T>G';'T>A'};
XY.alt_PRO=zeros(8,12);
XY.alt_CON=zeros(8,12);
XY.ref_PRO=zeros(8,12);
XY.ref_CON=zeros(8,12);


% maf
for b1=1:4, for b2=1:4, for s=1:6
            %keyboard
            c=cellstr([char(base(b1)) '-' char(snp(s)) '-' char(base(b2))]);
            k=find(ismember(X1.BRAB,c{1}));
            y=b2+4*mod(7-s-1,3);
            x=b1+4*floor(s/4);
            if (~isempty(k))
                XY.alt_PRO(x,y)=X1.PRO_ALT_BASES(k);
                XY.alt_CON(x,y)=X1.CON_ALT_BASES(k);
                XY.ref_PRO(x,y)=X1.PRO_REF_BASES(k);
                XY.ref_CON(x,y)=X1.CON_REF_BASES(k);
            end
            XY.col(x,y,:)=colors(s,:);
            XY.cat(x,y)=c;
            XY.snp(x,y)=snp(s);
        end,end,end

k=find(XY.ref_PRO<1)
XY.ref_PRO(k)=1
k=find(XY.ref_CON<1)
XY.ref_CON(k)=1




xp=0.2
    
if ismember(P.PROCON,'BOTH')
     z=(XY.alt_PRO+XY.alt_CON)./(XY.ref_PRO+XY.alt_PRO+XY.ref_CON+XY.alt_CONO)   
end
if ismember(P.PROCON,{'PRO'})
    z=XY.alt_PRO./(XY.ref_PRO+XY.alt_PRO)
end
if ismember(P.PROCON,{'CON'})
    z=XY.alt_CON./(XY.ref_PRO+XY.alt_PRO)
end

%z=100*z;

% PRO errors - relative to lego label convention (REF A OR C)
%zmax=1.1*max([XY.alt_PRO(:)./(XY.ref_PRO(:)+XY.alt_PRO(:)); XY.alt_CON(:)./(XY.ref_CON(:)+XY.alt_CON(:))])
zmax=1.15*max(z(:));
if length(P.zmax)>0
    zmax=P.zmax
end

%%
%clf;subplot(2,3,4)
h=bar3_with_colors(z,XY.col)
zlim([0 zmax])
%zlabel('error rate')
hz=zlabel('fraction of mismatched bases')
%hz=zlabel('percent of mismatched bases')
hz.FontSize=12;
set(gca,'xticklabel',[],'yticklabel',[],'ylim',[0.5 8.5]);
hb=gca;
%text(xp,1.5,['PRO REF>ALT'],'units','normalized','fontsize',11,'interpreter','none');
%text(xp,1.4,sprintf('n=%d',NS),'units','normalized','fontsize',11)
if P.baseline~=0
    for ih=1:length(h)
        z=get(h(ih),'ZData');
        z(z==0)=1e-9;
        set(h(ih),'ZData',z);
    end
end
if P.log
    set(gca,'zscale','log');
end
if length(P.ztick)>0
    set(gca,'ztick',P.ztick);
end

set(hb,'Position',get(hb,'Position')+[-0.01 0.0 0.0 0.02])

nc=zeros(6,1);
%e=XY.alt_PRO./(XY.ref_PRO+XY.alt_PRO);
for i=1:length(snp)
    k=strfindk(XY.snp(:),snp{i});
    nc(i)=sum(z(k))/sum(z(:));
end
nc=nc+1e-50;
hp=axes('Position', (get(hb,'Position')+[0.00,0.25,0, 0]).*[1 1 0.3 0.3]);
h=pie(nc,repmat({''},length(nc),1));
for i=1:6
    set(h(2*i-1),'facecolor',colors(i,:))    
end

h1=legend(hp,snp); %regexprep(snp,'>','>'));
set(h1,'linewidth',1,'FontSize',10)
%keyboard


%set(h1,'Position',get(hb,'Position')+[-0.075 -0.05 -0.05 0])
set(h1,'Position', (get(hb,'Position')+[0.2,0.2,0, 0]).*[1 1 0.04 0.3]);

%
%set(h1,'Position',[.42 .6 .07 .2],'linewidth',1,'FontSize',10)
q=regexprep(XY.cat(1:4,1:4),'-','_');
box1=cellfun(@(x) [x([1 2 7]) ' '],q,'UniformOutput',false);

hl=axes('Position', (get(hb,'Position')+[0.1,-0.11,0, 0]).*[1 1 0.75 0.95]);
h=bar3(ones(4,4),'w');
set(h,'EdgeAlpha',0.25);
set(gca,'zlim',[0 50]);
for i=1:4,for j=1:4
        hh(i,j)=text((j-0.35),(i)+0.15,4,box1(i,j),'horizontalAlign','left','interpreter','none'); %,'units','normalized')
    end,end
set(hl,'xtick',[],'ytick',[],'ztick',[],'visible','off');set(gcf,'color',[1 1 1]);




