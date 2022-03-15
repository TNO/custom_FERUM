function plot(a,d)
%
% Usage :  
%           plot(a)
%
% a -- The algorithm already trained
%      The colormap can be changed afterwards.
%
% Example:
%   d=gen(spiral({'n=1','m=50'}));
%   [r s0]=train(svm(kernel('rbf',1)),d)
%   plot(s0);  



%clf
%d=a.keep;
x=get_x(d);
y=get_y(d);
ax1=min(x); ax2=max(x);


granul=100;
minX=floor(ax1(1));
maxX=ceil(ax2(1));
minY=floor(ax1(2));
maxY=ceil(ax2(2));
axis_sz=[minX maxX minY maxY];
minX=minX-2; maxX=maxX+2; minY=minY-2; maxY=maxY+2;
gridx=[]; gridy=[];

mx=zeros(granul*granul,2);
for i=1:granul
 for j=1:granul
   mx((i-1)*granul+j,:)= [minX+(i-1)*(maxX-minX)/granul minY+(j-1)*(maxY-minY)/granul ] ;
   gridx(i)=minX+(i-1)*(maxX-minX)/granul;
   gridy(j)=minY+(j-1)*(maxY-minY)/granul;
end
end


temp=zeros(granul,granul);
dx=data(mx);
N=get_dim(dx);
blocks=N/50;
rr=[];
last=1;

% a.algorithm.use_signed_output=0;
%a.return_indices=1;
allx=[];
for i=blocks:blocks:N
       resX=test(a,get(dx,[last:i]));
       rr=[rr;resX.X'];allx=[allx;resX.X];
       last=i+1;
end

for i=1:granul
 for j=1:granul
  temp(i,j)= rr( (i-1)*granul+j);
  end
end
hold on;
%surf(1:granul,1:granul,temp'-1000);

% Ipos=find(allx(:,1)==+1);
% Ineg=find(allx(:,1)==-1);

% plot(mx(Ipos,1),mx(Ipos,2),'r.');
% plot(mx(Ineg,1),mx(Ineg,2),'b.');


%clf

% colormap('cool')
 pcolor(gridx, gridy, temp) ;
 shading interp;


%surf(gridx,gridy,temp');
%view(2);
%shading interp; 
%[c,h]=contour(gridx,gridy,temp',[0 0],'k');
%   clabel(c,h);
colorbar

if(1)

    Ipos=find(y==+1);
    Ineg=find(y==-1);
    h=plot(x(Ipos,1),x(Ipos,2),'rx'); hold on;
    set(h,'LineWidth',2,'MarkerSize',5);

    h=plot(x(Ineg,1),x(Ineg,2),'ko'); hold on;
    set(h,'LineWidth',2,'MarkerSize',5);
end

axis(axis_sz);






