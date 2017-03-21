% extract data from cellprofiler DefaultOUT.mat

% make sure the data file is in the current matlab folder
[NFkB]=importdata('DefaultOUTmat_Nuclei.xls');
% col 59 contains the cell id of the tracking
celldata=zeros(max(NFkB.data(:,59)),90);
t=0;
% reorder data to get time course of single cells
for findset=4:7804
    if length(NFkB.textdata{findset,1})>0
        t=t+1
    end
    celldata(NFkB.data(findset,59),t)=NFkB.data(findset,46);
end

%celldata=celldata(celldata>0)
celldata(celldata==0)=NaN

%% number of missing datapoints per cell
missnum=zeros(size(celldata,1),1)
for i=1:size(celldata,1)
    missnum(i)=sum(isnan(celldata(i,:)));
end
[B,I]=sort(missnum)

surf(celldata(I(1:80),1:80)')
plot(celldata(I(1:80),1:80)')

%% response integral
for i=1:size(celldata(I,:),1)
    datsum(i)=nansum((celldata(I(i),:)));
end
[M,N]=sort(datsum)
surf(celldata(N(end-80:end),1:80)')
%%
plot(celldata(N(end-80:end),1:80)')

%% response outliers
clear celldat sel
o=0
for i=1:size(celldata(I,:),1)
    sel(i)=sum(celldata(I(i),:)'>1.0e+04 *4)
    if sum(celldata(I(i),:)'>1.0e+04 *4)<1
        o=o+1
    celldat(i,1:80)=celldata(I(i),1:80);
    end
end
celldat(celldat==0)=NaN
figure; plot(celldat')
%% Heatmap
cgo=HeatMap(celldata(I([1:42 46:50]),1:79),'ColumnLabels',[0:78]*4);
set(cgo,'Colormap',jet);
addXLabel(cgo,'time in min','Fontsize',12)
addYLabel(cgo,'cellID','Fontsize',12)
addTitle(cgo,'NFkB nuclear localisation')
imwrite(cgo,'NFkBkymograph.tif')
%% Heatmap
cgo=HeatMap(celldat(1:50,:),'ColumnLabels',[0:79]*4);
set(cgo,'Colormap',jet);