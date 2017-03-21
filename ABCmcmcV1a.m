% ABC MCMC algorithm
close all
clear all

tic
itmax=50000; %number of iterations 20000=12min
D=zeros(1,itmax);
alpha=zeros(itmax,1);
Dat=cell(1,itmax);
ParDat=cell(1,itmax);
par=zeros(18,itmax);
par(1:18,1)=ones(18,1)*10^(-8);
%par(1:18,1)=1./(10.^(4*(rand(18,1)-0.5*ones(18,1)))*10^(8)); %initial parameters
rpar(1:18,1) = [5*10^(-9) 10^(-9) 10^(-8) 10^(-8) 5*10^(-9) 5*10^(-10) 10^(-7) 10^(-8) 10^(-7) 10^(-6) 10^(-8) 10^(-9) 10^(-8) 10^(-9) 10^(-8) 10^(-9) 10^(-8) 10^(-9)];

%par(1:18,1)=awgn(par(1:18,1),160);
k=zeros(18,itmax);
%k(1:18,1)=par(1:18,1);
var=10^(-8); % pertubation kernel variance
tmax=10^16; %ode solver integration time max
idxmax=65; %datapoints to 'regress'
        load('simdata'); % T0, Y0, idx
        sparT=T0(idx(1:idxmax));
        conDat=awgn(Y0(idx(1:idxmax),1)',1000);
            for m=2:8
                conDat=horzcat(conDat,abs(awgn(Y0(idx(1:idxmax),m)',1000)));
            end
%             D=1; i=0;
% while D>0.3
%     i=i+1;
for i=1:itmax
    D(i)=30;
    
    while D(i)>29
%         if i==1
%             par(1:18,1)=1./(10.^(4*(rand(18,1)-0.5*ones(18,1)))*10^(8)); %initial parameters
%         end
       k(1:18,i)=normrnd(par(1:18,i),par(1:18,i)/(10*log10(i)));%/(10*log10(i+1))); %draw candidate parameter
%          k(1:18,i)=awgn(par(1:18,i),180);
        % par(1:18,i)=i*ones(18,1);
        % k(1:18,i)=normrnd(par(1:18,i),5)*10^(-9);
        %x = [0:0.1:100];
        %figure;
        %y = normpdf(x,par(i),5);
        %plot(y)
        y0=zeros(8,1);
        [T0,Y0]=ode15s(@(t,y) Pi3kAktModelv8a_k(t,y,k(1,i),k(2,i),k(3,i),k(4,i),k(5,i),k(6,i),k(7,i),k(8,i),k(9,i),k(10,i),k(11,i),k(12,i),k(13,i),k(14,i),k(15,i),k(16,i),k(17,i),k(18,i)),[1 tmax],y0);
        Dat{i}.t=T0;
        Dat{i}.y=Y0;
        y0=zeros(8,1);
        [tt0,yy0]=ode15s(@(t,y) Pi3kAktModelv8a_k(t,y,par(1,i),par(2,i),par(3,i),par(4,i),par(5,i),par(6,i),par(7,i),par(8,i),par(9,i),par(10,i),par(11,i),par(12,i),par(13,i),par(14,i),par(15,i),par(16,i),par(17,i),par(18,i)),[1 tmax],y0);
        ParDat{i}.t=tt0;
        ParDat{i}.y=yy0;

        
        %idx ..index of time data closest to equaly spaced sparse data time vektor
        %concatenate sparse concentration time traces of all species

        load('simdata'); % T0, Y0, idx

        a=vdsearchn(Dat{i}.t',sparT); %find index 'a' closest to sparT
        concata=Dat{i}.y(a(1:idxmax),1)';
            for n=2:8
                %conDat.t=horzcat(conDat.t,Dat(j).t')
                concata=horzcat(concata,Dat{i}.y(a(1:idxmax),n)');
            end
        conSima(i,:)=concata;
        
        b=vdsearchn(ParDat{i}.t',sparT); %find index 'a' closest to sparT
        concatb=ParDat{i}.y(b(1:idxmax),1)';
            for n=2:8
                %conDat.t=horzcat(conDat.t,Dat(j).t')
                concatb=horzcat(concatb,ParDat{i}.y(b(1:idxmax),n)');
            end
        conSimb(i,:)=concatb;
        


        [ID,D(i)] =knnsearch(conSima,conDat);
         [D(i), i/10^4]
    end
    %plot(conSim');
%     plot(conSim(ID,:)); hold on;
%     plot(conDat); hold off;
    x=1;
    for l=1:idxmax
        x = x* mvnpdf(Y0(idx(l),:),Dat{i}.y(a(l),:),eye(8,8)*10^(-log2(i)*0.1));%*10^(-round(i/3000))); %p(data|k-sim)
    end
    y=1;
    for l=1:idxmax
        y = y* mvnpdf(Y0(idx(l),:),ParDat{i}.y(b(l),:),eye(8,8)*10^(-log2(i)*0.1));%*10^(-round(i/3000))); %p(data|par-sim)
    end
    alpha(i)=min(1,y/x);
    if alpha(i)<1 %set par=k with probability alpha
        if alpha(i)>rand
        par(1:18,i+1)=k(1:18,i);
%         disp('parchange')
        else
        par(1:18,i+1)=par(1:18,i);
        end
    else
        par(1:18,i+1)=par(1:18,i);
    end
% x
% y
end
toc
% u=figure;plot(par(:,1:i)'); hold on; title('Accepted Parameters')
% imwrite(u)
v=figure;plot(k(:,1:i)'); hold on; title('Proposed Parameters')
savefig(v,'SC170320ProposedParameters')
w=figure;
    plot(conSima(i,:)); hold on; title('Merged Time-course Vectors Data vs. Simulation ')
    plot(conSima(1,:),'o');
%     plot(conSima(1,:),'x');
    plot(conDat,'x'); hold off;
savefig(w,'SC170320timecourses')

z=3; %zeilen
s=round(size(par,1)/z); %spalten
his=figure;
for p=1:size(par,1)
subplot(z,s,p);
hist(k(p,round(itmax*0.9):round(itmax)),50); hold on; title(num2str(rpar(p,1)));
end
savefig(his,'SC170320histograms')

dpl=figure; plot(D(2:end));
savefig(dpl,'SC170320convergence')

save('SC170320abcmcmcAdaptiveL','D','par','k','ParDat','Dat','conDat','conSima','conSimb')