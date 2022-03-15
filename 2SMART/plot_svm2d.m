function plot_svm2d(SVC,gca_param,varargin)

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
   ax1 = [ gca_param.XLim(1) gca_param.YLim(1) ]; ax2 = [ gca_param.XLim(2) gca_param.YLim(2) ];
end

minX = floor(ax1(1));
maxX = ceil(ax2(1));
minY = floor(ax1(2));
maxY = ceil(ax2(2));


granul = 300+1;

[X,Y] = meshgrid(linspace(minX,maxX,granul),linspace(minY,maxY,granul));
mx = [ X(:) Y(:) ];

sflag = SVC.algorithm.use_signed_output;
SVC.algorithm.use_signed_output = 0;

resX = test(SVC,data(mx));
resX = sign(get_x(resX)) .* sqrt(abs(get_x(resX)));

C = reshape(resX,granul,granul);

if nargin == 13 | nargin == 7
   [ G, dummy ] = gfunbasic(1,u_to_x(mx',probdata),'no ',probdata,analysisopt,gfundata);
   Cref = reshape(G,granul,granul);   
end

hold on;
colormap('gray')
%map = colormap('jet');
%colormap(flipud(map));
pcolor(X,Y,C) ;
%colorbar
shading interp

if nargin == 8 | nargin == 13

   if ~isempty(Usvm)
      plot(Usvm(1,:),Usvm(2,:),'g.','Markersize',4)
   end

   if n_pt_supp_clust > 0
      plot(Usvm_in_margin(1,:),Usvm_in_margin(2,:),'y.','Markersize',4)
   else
      plot(Usvm_in_margin_sorted(1,:),Usvm_in_margin_sorted(2,:),'y.','Markersize',4)
   end
   if n_pt_supp_switch > 0
      plot(Usvm_switch(1,:),Usvm_switch(2,:),'m.','Markersize',4)
   end

end

[c,hs] = contour(X,Y,C,[0 0],'k--','LineWidth',2);
[c,hr] = contour(X,Y,C,[-1 -1],'r--','LineWidth',2);
[c,hb] = contour(X,Y,C,[1 1],'b--','LineWidth',2);
if nargin == 13 | nargin == 7
   for istep = 0:(Nb_step-1)
      [c,hgth] = contour(X,Y,Cref,[G0(istep+1) G0(istep+1)],'k--','LineWidth',2);
      hgth(istep+1) = hgth;
   end
   [c,hgth] = contour(X,Y,Cref,[G0(Nb_step+1) G0(Nb_step+1)],'r-','LineWidth',2);
   hgth(Nb_step+1) = hgth;
   [c,hg0] = contour(X,Y,Cref,[0 0],'k-','LineWidth',2);
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
   h = plot(x(i,1),x(i,2),col);
   alphasv = ceil((abs(SVC.alpha(i))/amax)*4);
   set(h,'LineWidth',alphasv,'MarkerSize',8+alphasv);
end


axis([ minX maxX minY maxY ]);