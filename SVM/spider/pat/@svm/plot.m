function data_square_plot(a,d)
%
% Plots a surface plot of an svm which uses !_2_! dimensional inputs.
% Red decision line is drawn.
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
if(nargin <2)
    d=a.Xsv;
end

x=d.X;
y=d.Y;
ax1=min(x); ax2=max(x);
ax1=[-5 -5], ax2=[5 5]
expand = 0.2 * [ax2-ax1];
ax1 = ax1 - expand;
ax2 = ax2 + expand;

c1=find(y==1);
c2=find(y==-1);



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

sflag=a.algorithm.use_signed_output;
a.algorithm.use_signed_output=0;

resX=test(a,data(mx));
resX=sign(get_x(resX)).*sqrt(abs(get_x(resX)));

for i=1:granul
    for j=1:granul
        temp(i,j)= resX( (i-1)*granul+j);
    end
end

hold on;

% FeatureLines = [0 -1 1]';  %cheap hack to only get the decision boundary
colormap('gray')
%map = colormap('jet');
%colormap(flipud(map));
pcolor(gridx, gridy, temp') ;
shading interp
[c,hs] = contour(gridx, gridy, temp',[0 0],'k--','LineWidth',2);
[c,hr] = contour(gridx, gridy, temp',[-1 -1],'r--','LineWidth',2);
[c,hb] = contour(gridx, gridy, temp',[1 1],'b--','LineWidth',2);

%colorbar


sv=find(abs(a.alpha)>1e-7);


%sv=find(abs(a.alpha)>max(abs(a.alpha))/100);
amax=max(abs(a.alpha));

[r]=test(a,a.Xsv);

if(nargin <2)

    if 1
        for i=sv'
            col='co';
            if (r.X(i,:)*r.Y(i,:)<0.9) col='cs';   end;  %% margin error
            h=plot(x(i,1),x(i,2),col);
            alpha=ceil((abs(a.alpha(i))/amax)*4);
            set(h,'LineWidth',alpha,'MarkerSize',8+alpha);
        end
    end

else
    
    Ip=find(y==+1);
    In=find(y==-1);
    
     plot(x(Ip,1),x(Ip,2),'r.','MarkerSize',2)
     plot(x(In,1),x(In,2),'b.','MarkerSize',2)
    
    Z=get_x(a.Xsv);
    plot(Z(:,1),Z(:,2),'go','MarkerSize',2)
end

axis(axis_sz);

