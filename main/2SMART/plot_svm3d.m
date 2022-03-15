function plot_svm3d(SVC,gca_param,varargin)

if nargin == 8
   Usvm                  = varargin{1};
   n_pt_supp_clust       = varargin{2};
   Usvm_in_margin        = varargin{3};
   Usvm_in_margin_sorted = varargin{4};
   n_pt_supp_switch      = varargin{5};
   Usvm_switch           = varargin{6};
elseif nargin == 13
   Usvm                  = varargin{1};
   n_pt_supp_clust       = varargin{2};
   Usvm_in_margin        = varargin{3};
   Usvm_in_margin_sorted = varargin{4};
   n_pt_supp_switch      = varargin{5};
   Usvm_switch           = varargin{6};
   probdata              = varargin{7};
   analysisopt           = varargin{8};
   gfundata              = varargin{9};
   Nb_step               = varargin{10};
   G0                    = varargin{11};
elseif nargin == 7
   probdata              = varargin{1};
   analysisopt           = varargin{2};
   gfundata              = varargin{3};
   Nb_step               = varargin{4};
   G0                    = varargin{5};
end

D = SVC.Xsv;

x = D.X;
y = D.Y;

if gca_param.auto_gca
   ax1 = min(x); ax2 = max(x);
else
   ax1 = [ gca_param.XLim(1) gca_param.YLim(1) gca_param.ZLim(1) ]; ax2 = [ gca_param.XLim(2) gca_param.YLim(2) gca_param.ZLim(2) ];
end

minX = floor(ax1(1));
maxX = ceil(ax2(1));
minY = floor(ax1(2));
maxY = ceil(ax2(2));
minZ = floor(ax1(3));
maxZ = ceil(ax2(3));


granul = 100+1;

[X,Y,Z] = meshgrid(linspace(minX,maxX,granul),linspace(minY,maxY,granul),linspace(minZ,maxZ,granul));
mx = [ X(:) Y(:) Z(:) ];

sflag = SVC.algorithm.use_signed_output;
SVC.algorithm.use_signed_output = 0;

resX = test(SVC,data(mx));
resX = sign(get_x(resX)) .* sqrt(abs(get_x(resX)));

C = reshape(resX,granul,granul,granul);

if nargin == 13 | nargin == 7
   [ G, dummy ] = gfunbasic(1,u_to_x(mx',probdata),'no ',probdata,analysisopt,gfundata);
   Cref = reshape(G,granul,granul,granul);   
end

hold on;

P = patch(isosurface(X,Y,Z,C,0));
isonormals(X,Y,Z,C,P);
set(P,'facecolor',[0.5 0.5 0.5],'edgecolor','none');

daspect([1 1 1]);
view(3); axis tight; grid on;

if nargin == 8 | nargin == 13

%      plot(Usvm(1,:),Usvm(2,:),Usvm(3,:),'g.','Markersize',4)

      if n_pt_supp_clust > 0
         plot3(Usvm_in_margin(1,:),Usvm_in_margin(2,:),Usvm_in_margin(3,:),'y.','Markersize',4)
      else
         plot3(Usvm_in_margin_sorted(1,:),Usvm_in_margin_sorted(2,:),Usvm_in_margin_sorted(3,:),'y.','Markersize',4)
      end
      if n_pt_supp_switch > 0
         plot3(Usvm_switch(1,:),Usvm_switch(2,:),Usvm_switch(3,:),'m.','Markersize',4)
      end

end

if nargin == 13 | nargin == 7
   Pgth = patch(isosurface(X,Y,Z,Cref,G0(Nb_step+1)));
   isonormals(X,Y,Z,Cref,Pgth);
   set(Pgth,'facecolor',[1 0 0],'edgecolor','none');
end

sv = find(abs(SVC.alpha)>1e-7);
%sv = find(abs(SVC.alpha)>max(abs(SVC.alpha))/100);
amax = max(abs(SVC.alpha));
r = test(SVC,SVC.Xsv);
for i = sv'
   col = 'co';
   if r.X(i,:)*r.Y(i,:) < 0.9
      col = 'cs';
   end;  %% margin error
%    h = plot3(x(i,1),x(i,2),x(i,3),col);
%    alpha = ceil((abs(SVC.alpha(i))/amax)*4);
%    set(h,'LineWidth',alpha,'MarkerSize',8+alpha);
   alphasv = ceil((abs(SVC.alpha(i))/amax)*4);
   [Xsv,Ysv,Zsv] = sphere(20);
   Xsv = x(i,1)+(6+alphasv)/50*Xsv;
   Ysv = x(i,2)+(6+alphasv)/50*Ysv;
   Zsv = x(i,3)+(6+alphasv)/50*Zsv;
   hsv = surf(Xsv,Ysv,Zsv);
   set(hsv,'FaceColor','c','EdgeColor','none');
end

camlight; 
lighting gouraud;

alpha(P,.9)
if nargin == 13 | nargin == 7
   alpha(Pgth,.6)
end

%axis([ minX maxX minY maxY minZ maxZ ]);

