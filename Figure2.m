
clear
rmpath /Users/stewart/CancerGenomeAnalysis/trunk/matlab/mike
cd ~/GoogleDrive/Comment_to_Science_damage/data/
X1=load('~/GoogleDrive/Comment_to_Science_damage/data/TCGA.1933.error_rate.29Oct2017.mat')
%printStruct(X1,-1,'~/GoogleDrive/Comment_to_Science_damage/data/TCGA.1933.error_rate.30Oct2017.txt')

X1.error_rate=X1.total_alt_bases./X1.total_bases;
[q,k]=min(X1.error_rate)
look(X1,k)



subplot('position',[0.1 0.1 0.23 0.4])
plot(X1.GIV,X1.oxoG_alt_bases./X1.total_bases,'rx')
hold on
plot(X1.GIV,X1.nonoxoG_alt_bases./X1.total_bases,'b+')
plot(X1.GIV,X1.oxoG_alt_bases./X1.total_bases,'rx')
hold off
xlim([0 35])
xlabel('GIV_G_T corrected score')
xlabel({'GIV_{G\_T} corrected score'},'interpreter','tex')
ylabel('Errors per bases sequenced')
legend({'oxoG errors','non-oxoG errors'},'location','NW')
grid on
yt=0:0.0005:0.0025; set(gca,'ytick',yt,'yticklabels',yt);ylim([0 3e-3])
text(-0.29,1.01,'C','Units', 'Normalized', 'VerticalAlignment', 'Top','fontsize',18)

subplot('position',[0.1 0.79 0.23 0.14])
xb=0:35
nb=hist(X1.GIV,xb)
h=bar(xb,nb,'hist')
xlim([0 35])
%set(gca,'yscale','log')
set(gca,'xticklabel',[])
set(h,'FaceColor',0.5*[1 1 1])
h.Vertices(h.Vertices==0)=0.9
set(gca,'ytick',200:200:800)
ylim([0 800])
ylabel('Tumors')
grid on
text(-0.29,1.1,'A','Units', 'Normalized', 'VerticalAlignment', 'Top','fontsize',18)

subplot('position',[0.1 0.52 0.23 0.25])
plot(X1.GIV,100*X1.oxoG_alt_bases./X1.total_alt_bases,'k.')
xlim([0 35])
ylabel({'Percentage of errors ','from 8-oxoG damage'})
grid on
ylim([0 100])
set(gca,'xticklabel',[])
yt=20:20:100; set(gca,'ytick',yt,'yticklabels',yt);
text(-0.29,1.1,'B','Units', 'Normalized', 'VerticalAlignment', 'Top','fontsize',18)
set(gcf,'Position',[100 200 1100 600])

% D
cd /Users/stewart/GoogleDrive/Comment_to_Science_damage/data
I = imread('/Users/stewart/GoogleDrive/Comment_to_Science_damage/data/Chen.Fig4E.png');

subplot(2,3,2)
%imshow(I)
p=double(I(:,:,1)-I(:,:,2)-I(:,:,3))./sum(I,3)

z1=flipud(1-(p>0.05));
%imagesc(z1)
d=size(z1)
s0=[1 d(1) 1 d(2)]
%s1=[0 100 1 7.8]
y0=(1:d(1))
y1=(y0-1)*(100/(d(1)-1))
x0=(1:d(2))
x1=(x0-1)*(6.65/(d(2)-1))+1

z=flipud(I);
imagesc(x1,y1,z)
axis xy
h=xlabel('GIV_{G\_T}'); set(h,'interpreter','tex')
ylabel({'percentage of estimated',' false positive somatic variants'})
%set(gcf,'position',[50 50 1000 1200])

%text(0.02,0.98,'A','Units', 'Normalized', 'VerticalAlignment', 'Top','fontsize',18)
text(-0.29,1.05,'D','Units', 'Normalized', 'VerticalAlignment', 'Top','fontsize',18)

% E
P=load_tsv('/Users/stewart/GoogleDrive/Cancer/damage/mc3/mc3.oxoG.damage.GIV.FDR.prior.pairs.txt')

SET=flipud(unique(P.SET)); n=length(SET)

v='includePAAD'
if (1)
    v='noPAAD'
    P=trimStruct(P,~ismember(P.SET,'PAAD'))
    SET=flipud(unique(P.SET)); n=length(SET)
end
SET{1}='THCA'
SET{2}='UCS'

subplot(2,3,3)

plot(P.GIV_GT_QUALITYSCORE20,100*P.fdr_part_posterior,'ro','markersize',4,'markerfacecolor','r')
grid on
xlim([1 7.5])
ylim([0 100])
ylabel('MuTect FDR (%) BEFORE filter .')
%xlabel({'GIV_G_T corrected'},'interpreter','none')
xlabel({'GIV_{G\_T} corrected'},'interpreter','tex')
x = [1 7.5 7.5 1 1];
y = [0 0 100 100 0];
p=patch(x,y,0.7*[1 1 1] );
set(p,'FaceAlpha',0.5);
set(p,'EdgeColor','none');
%r = rectangle('Position',[1 0 7.5 100]')
%r.EdgeColor=[0 0 0]
%r.FaceColor=0.8*[1 1 1]
hold on;
plot(P.GIV_GT_QUALITYSCORE20,100*P.fdr_part_posterior,'ro','markersize',4,'markerfacecolor','r')
hold off
grid on;
text(-0.29,1.05,'E','Units', 'Normalized', 'VerticalAlignment', 'Top','fontsize',18)

% F
subplot(2,3,5)

plot(P.GIV_GT_QUALITYSCORE20,100*P.fdr_part_posterior,'.','markersize',11)
grid on
xlim([0 35])
for i=2:n
    hold on
    k=ismember(P.SET,SET(i));
    plot(P.GIV_GT_QUALITYSCORE20(k),100*P.fdr_part_posterior(k),'.','markersize',11)
end
hold off
ylabel('MuTect FDR (%) BEFORE filter ')
%xlabel({'GIV_G_T  corrected'},'interpreter','none')
xlabel({'GIV_{G\_T} corrected'},'interpreter','tex')
grid on;
[h,icons,plots,legend_text]=legend(SET,'location','SE')
NSETS=length(SET)
kl=findobj(icons,'type','Line')
for o=kl(2:2:(2*NSETS))'
    o.MarkerSize=20
end
kt=findobj(icons,'type','Text')
for o=kt'
    o.FontSize=14
end
text(-0.29,1.1,'F','Units', 'Normalized', 'VerticalAlignment', 'Top','fontsize',18)

%G
subplot(2,3,6)
plot(P.GIV_GT_QUALITYSCORE20,100*P.fdr_part_posterior_pass,'.','markersize',11)
grid on
xlim([0 35])
for i=2:n
    hold on
    k=ismember(P.SET,SET(i));
    plot(P.GIV_GT_QUALITYSCORE20(k),100*P.fdr_part_posterior_pass(k),'.','markersize',11)
end
hold off
ylabel('MuTect FDR (%) AFTER filter')
xlabel({'GIV_{G\_T} corrected'},'interpreter','tex')
grid on;

[h,icons,plots,legend_text]=legend(SET,'location','SE')
NSETS=length(SET)
kl=findobj(icons,'type','Line')
for o=kl(2:2:(2*NSETS))'
    o.MarkerSize=20
end
kt=findobj(icons,'type','Text')
for o=kt'
    o.FontSize=14
end
ylim([0 100])
text(-0.29,1.1,'G','Units', 'Normalized', 'VerticalAlignment', 'Top','fontsize',18)

% plot
fplt=['~/GoogleDrive/Comment_to_Science_damage/matlab/plots/Fig2.' TODAY ]
saveas(gcf,[fplt '.png'],'png')
print([fplt '.png'],'-dpng')
print([fplt '.eps'],'-depsc','-painters')
print([fplt '.svg'],'-dsvg','-painters')

%saveas(gcf,[fplt '.eps'],'epsc')
%saveas(gcf,[fplt '.svg'],'svg')


