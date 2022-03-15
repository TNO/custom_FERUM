function data = subset_subset_graph(data,option)

plot_small = 1;
if plot_small
   blueMarker = '.';   blueMarkersize = 10; blueLineWidth = 2; 
   redMarker  = '^';   redMarkersize  = 4;  redLineWidth  = 1; 
   genMarker  = '*';   genMarkersize  = 4;  genLineWidth  = 0.5; 
else
   blueMarker = 'x';   blueMarkersize = 8;  blueLineWidth = 1.5; 
   redMarker  = '^';   redMarkersize  = 8;  redLineWidth  = 2.5; 
   genMarker  = '*';   genMarkersize  = 8;  genLineWidth  = 1.5; 
end

auto_gca = 0;
if ~auto_gca
   XLim = [-5 5];
   YLim = [-5 5];
   XTick = [-5 0 5];
   YTick = [-5 0 5];
end

% plot_thresholds = 0;
% if plot_thresholds
%     "..."
%     u1gth = linspace(XLim(1),XLim(2),1000);
%     for istep = 1:length(data.y)
%        u2gth(istep,:) = "g(u1gth)" - data.y(istep);
%     end
%     u1g0 = u1gth;
%     u2g0 = "g(u1g0)";
% end

switch option
    
   case 0 % Samples lower than threshold value
       
      figid = 2*data.Nb_step+1;
       
      if data.Nb_step == 0
         data.U1Lim = [ 0 0 ];
         data.U2Lim = [ 0 0 ];
      end
      
      U1Lim = [ min(data.U(1,data.Indices(end,1):data.Indices(end,2))) max(data.U(1,data.Indices(end,1):data.Indices(end,2))) ];
      U2Lim = [ min(data.U(2,data.Indices(end,1):data.Indices(end,2))) max(data.U(2,data.Indices(end,1):data.Indices(end,2))) ];
      if U1Lim(1) < data.U1Lim(1), data.U1Lim(1) = floor(U1Lim(1)*2)/2; end  
      if U1Lim(2) > data.U1Lim(2), data.U1Lim(2) = ceil(U1Lim(2)*2)/2; end
      if U2Lim(1) < data.U2Lim(1), data.U2Lim(1) = floor(U2Lim(1)*2)/2; end
      if U2Lim(2) > data.U2Lim(2), data.U2Lim(2) = ceil(U2Lim(2)*2)/2; end
      for i = 1:figid-1
         if auto_gca
            set(gca(i),'XLim',data.U1Lim,'YLim',data.U2Lim);
         end
         axis square
      end
      
      figure(figid)
      title([ 'Step #' num2str(data.Nb_step) ],'FontSize',14)
      hold on
      h1 = plot(data.U(1,data.Indices(end,1):data.Indices(end,2)),data.U(2,data.Indices(end,1):data.Indices(end,2)),blueMarker);
      set(h1,'Markersize',blueMarkersize,'LineWidth',blueLineWidth,'Color','b','MarkerFaceColor','b');
      imax_Indgerm = find(data.Indgerm(end,:)==0);
      if isempty(imax_Indgerm)
         imax_Indgerm = size(data.Indgerm(end,:),2);
      else
         imax_Indgerm = imax_Indgerm(1)-1;
      end
      h2 = plot(data.U(1,data.Indgerm(end,1:imax_Indgerm)),data.U(2,data.Indgerm(end,1:imax_Indgerm)),redMarker);
      set(h2,'Markersize',redMarkersize,'LineWidth',redLineWidth,'Color','r','MarkerFaceColor','r');
%       if plot_thresholds
%          hgth = plot(u1gth,u2gth(end,:),'-');
%          set(hgth,'LineWidth',2,'Color','r')
%          hg0 = plot(u1g0,u2g0,'-');
%          set(hg0,'LineWidth',2,'Color','k')
%       end
      if auto_gca
         set(gca,'FontSize',14,'XLim',data.U1Lim,'YLim',data.U2Lim);
      else
         set(gca,'FontSize',14,'XLim',XLim,'YLim',YLim,'XTick',XTick,'YTick',YTick);
         grid on
         box on
      end
      axis square
      
      xlabel('{\itu}_1','FontSize',14);
      ylabel('{\itu}_2','FontSize',14);
      
      disp('Press a key to continue')
      pause
      
    case 1 % MCMC, all states

      figid = 2*data.Nb_step;

      figure(figid)
      clf
      title([ 'Step #' num2str(data.Nb_step) ],'FontSize',14)
      hold on
      h1 = plot(data.U(1,data.Indices(end-1,1):data.Indices(end-1,2)),data.U(2,data.Indices(end-1,1):data.Indices(end-1,2)),blueMarker);
      set(h1,'Markersize',blueMarkersize,'LineWidth',blueLineWidth,'Color','b','MarkerFaceColor','b');
      h3 = plot(data.U(1,data.Indices(end,1):data.Indices(end,2)),data.U(2,data.Indices(end,1):data.Indices(end,2)),genMarker);
      set(h3,'Markersize',genMarkersize,'LineWidth',genLineWidth,'Color','g');
      imax_Indgerm = find(data.Indgerm(end,:)==0);
      if isempty(imax_Indgerm)
         imax_Indgerm = size(data.Indgerm(end,:),2);
      else
         imax_Indgerm = imax_Indgerm(1)-1;
      end
      h2 = plot(data.U(1,data.Indgerm(end,1:imax_Indgerm)),data.U(2,data.Indgerm(end,1:imax_Indgerm)),redMarker);
      set(h2,'Markersize',redMarkersize,'LineWidth',redLineWidth,'Color','r','MarkerFaceColor','r');
%       if plot_thresholds
%          hgth = plot(u1gth,u2gth(end,:),'-');
%          set(hgth,'LineWidth',2,'Color','r')
%          hg0 = plot(u1g0,u2g0,'-');
%          set(hg0,'LineWidth',2,'Color','k')
%       end
      if auto_gca
         set(gca,'FontSize',14,'XLim',data.U1Lim,'YLim',data.U2Lim);
      else
         set(gca,'FontSize',14,'XLim',XLim,'YLim',YLim,'XTick',XTick,'YTick',YTick);
         grid on
         box on
      end
      axis square
      
      xlabel('{\itu}_1','FontSize',14);
      ylabel('{\itu}_2','FontSize',14);
      
      disp('Press a key to continue')
      pause

   case 2 % MCMC, one single state
      
      figid = 2*data.Nb_step;

      figure(figid)
      hold on
      
      U1Lim = [ min(data.Usubtemp(1,:)) max(data.Usubtemp(1,:)) ];
      U2Lim = [ min(data.Usubtemp(2,:)) max(data.Usubtemp(2,:)) ];
      if U1Lim(1) < data.U1Lim(1), data.U1Lim(1) = floor(U1Lim(1)*2)/2; end  
      if U1Lim(2) > data.U1Lim(2), data.U1Lim(2) = ceil(U1Lim(2)*2)/2; end
      if U2Lim(1) < data.U2Lim(1), data.U2Lim(1) = floor(U2Lim(1)*2)/2; end
      if U2Lim(2) > data.U2Lim(2), data.U2Lim(2) = ceil(U2Lim(2)*2)/2; end
      for i = 1:figid-1
         if auto_gca
             set(gca(i),'XLim',data.U1Lim,'YLim',data.U2Lim);
         end
         axis square
         xlabel('{\itu}_1','FontSize',14);
         ylabel('{\itu}_2','FontSize',14);
      end
      
      if data.Nb_generation == 1
          
         title([ 'Step #' num2str(data.Nb_step)-1 ],'FontSize',14)
         h1 = plot(data.U(1,data.Indices(end,1):data.Indices(end,2)),data.U(2,data.Indices(end,1):data.Indices(end,2)),blueMarker);
         set(h1,'Markersize',blueMarkersize,'LineWidth',blueLineWidth,'Color','b','MarkerFaceColor','b');
         imax_Indgerm = find(data.Indgerm(end,:)==0);
         if isempty(imax_Indgerm)
            imax_Indgerm = size(data.Indgerm(end,:),2);
         else
            imax_Indgerm = imax_Indgerm(1)-1;
         end
         h2 = plot(data.U(1,data.Indgerm(end,1:imax_Indgerm)),data.U(2,data.Indgerm(end,1:imax_Indgerm)),redMarker);
         set(h2,'Markersize',redMarkersize,'LineWidth',redLineWidth,'Color','r','MarkerFaceColor','r');
         if auto_gca
            set(gca,'FontSize',14,'XLim',data.U1Lim,'YLim',data.U2Lim);
         else
            set(gca,'FontSize',14,'XLim',XLim,'YLim',YLim,'XTick',XTick,'YTick',YTick);
            grid on
            box on
         end
         axis square
         
         pause(1)
         
      end
      
      for i = 1:(data.Nb_generation-1)
         ht = data.ht;
         hthl = data.hthl;
         set(ht(i),'Markersize',genMarkersize,'LineWidth',genLineWidth,'Color','g');
%asl         if ~isnan(hthl(i))
         if ~isstruct(hthl(i))             
            set(hthl(i),'Markersize',genMarkersize,'LineWidth',genLineWidth,'Color','k');
         end
      end

      title([ 'Step #' num2str(data.Nb_step) ' - Generation #' num2str(data.Nb_generation) ])

      ht(data.Nb_generation) = plot(data.Usubtemp(1,:),data.Usubtemp(2,:),genMarker);
      set(ht(data.Nb_generation),'Markersize',genMarkersize,'LineWidth',genLineWidth,'Color','y');
      if ~isempty(data.ind)
         hthl(data.Nb_generation) = plot(data.Usubtemp(1,data.ind),data.Usubtemp(2,data.ind),genMarker);
         set(hthl(data.Nb_generation),'Markersize',genMarkersize,'LineWidth',genLineWidth,'Color','m');
      else
         hthl(data.Nb_generation) = nan;
      end
%       if plot_thresholds
%          hgth = plot(u1gth,u2gth(end,:),'-');
%          set(hgth,'LineWidth',2,'Color','r')
%          hg0 = plot(u1g0,u2g0,'-');
%          set(hg0,'LineWidth',2,'Color','k')
%       end
      if auto_gca
         set(gca,'FontSize',14,'XLim',data.U1Lim,'YLim',data.U2Lim);
      else
         set(gca,'FontSize',14,'XLim',XLim,'YLim',YLim,'XTick',XTick,'YTick',YTick);
         grid on
         box on
      end
      axis square
      
      data.ht = ht;
      data.hthl = hthl;
      
      xlabel('{\itu}_1','FontSize',14);
      ylabel('{\itu}_2','FontSize',14);
      
      pause(1)

    case 3 % Final plot (all steps)

      figid = 2*data.Nb_step+2;
      
      figure(figid)
      title('Final','FontSize',14)
      hold on

      map = colormap(jet(data.Nb_step+1));
      lestr = '{';
      for i = 0:data.Nb_step
         imax_Indgerm = find(data.Indgerm(i+1,:)==0);
         if isempty(imax_Indgerm)
            imax_Indgerm = size(data.Indgerm(i+1,:),2);
         else
            imax_Indgerm = imax_Indgerm(1)-1;
         end
         h(i+1) = plot(data.U(1,data.Indgerm(i+1,1:imax_Indgerm)),data.U(2,data.Indgerm(i+1,1:imax_Indgerm)),redMarker);
         set(h(i+1),'Markersize',redMarkersize,'LineWidth',redLineWidth,'Color',map(i+1,:),'MarkerFaceColor',map(i+1,:));
         lestr = [ lestr ' ''Step #' num2str(i) ''' ' ];
      end
      lestr = [ lestr '}' ];
%       if plot_thresholds
%          for istep = 1:length(data.y)
%             hgth(istep) = plot(u1gth,u2gth(istep,:),'-');
%             set(hgth(istep),'LineWidth',2,'Color',map(istep,:))
%          end
%       end
      if auto_gca
         set(gca,'FontSize',14,'XLim',data.U1Lim,'YLim',data.U2Lim);
      else
         set(gca,'FontSize',14,'XLim',XLim,'YLim',YLim,'XTick',XTick,'YTick',YTick);
         grid on
         box on
      end
      axis square
      eval([ 'legend(' lestr ');' ]);
      
      xlabel('{\itu}_1','FontSize',14);
      ylabel('{\itu}_2','FontSize',14);

end