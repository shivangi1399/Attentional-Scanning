function selected_saccades = selectSaccades_visual(saccades,n)
% SELECTSACCADES_VISUAL(saccades,n) Select n saccades from input saccades
    fprintf('\n[select_saccades_visual]\n')
    try
        assert(all(arrayfun(@(x) ~isempty(x.parent), saccades)));       
        assert(numel(unique(arrayfun(@(x) x.parent, saccades)))==1);
        plot_option = 1;
    catch
        warning('saccades have no common parent asigned')
        plot_option = 2; % WAIT FOR IMPLEMENTATION
    end
    
    switch plot_option
        case 1 % have common parent GazeEpoch
            pGazeEpoch = saccades(1).parent;           
            UserData.saccades   = saccades;
            UserData.activeInd  = 1;
            UserData.flagSelected    = zeros(size(saccades)); 
            UserData.flagProcessed = zeros(size(saccades));
            
            h.fig = figure('Name',sprintf('saccade %03d/%03d, %03d processed, %03d selected',...
                UserData.activeInd,numel(UserData.saccades), sum(UserData.flagProcessed), sum(UserData.flagSelected)),'OuterPosition',[500,300,800,600]);
            set(h.fig,'UserData',UserData);
            
            
            %init, x trace
            subplot(2,4,[1,2])
            ylim(pGazeEpoch.scene.limit(1:2));xlim([min(pGazeEpoch.time),max(pGazeEpoch.time)]);
            plot(pGazeEpoch.time,pGazeEpoch.position(:,1),'-','Color',[0.25,0.25,0.25],'LineWidth',0.5);hold on;
            plot([0 0],get(gca,'Ylim'),'k--');
            h.xtrace = plot(UserData.saccades(UserData.activeInd).time, UserData.saccades(UserData.activeInd).position(:,1),'-r','LineWidth',1);
            ylabel('Position');title('Horizontal component');grid minor;
            
            % init, y trace
            subplot(2,4,[5,6])
            ylim(pGazeEpoch.scene.limit(3:4));xlim([min(pGazeEpoch.time),max(pGazeEpoch.time)]);
            plot(pGazeEpoch.time,pGazeEpoch.position(:,2),'-','Color',[0.25,0.25,0.25],'LineWidth',0.5);hold on;
            plot([0 0],get(gca,'Ylim'),'k--');
            h.ytrace = plot(UserData.saccades(UserData.activeInd).time, UserData.saccades(UserData.activeInd).position(:,2),'-r','LineWidth',1);
            ylabel('Position');title('Vertical component');xlabel('Time [sec]');grid minor;
            
            % init, 2d trajectory
            subplot(2,4,[3,4,7,8])
            image(pGazeEpoch.scene.limit(1:2),pGazeEpoch.scene.limit(4:-1:3),...
                pGazeEpoch.scene.display);hold on;
            set(gca,'YDir','normal');
            set(gca,'Visible','off');
            plot(pGazeEpoch.position(:,1),pGazeEpoch.position(:,2),'-','Color',[0.25,0.25,0.25],'LineWidth',0.5);
            h.trajectory = plot(UserData.saccades(UserData.activeInd).position(:,1),UserData.saccades(UserData.activeInd).position(:,2),'.-r','LineWidth',1);
            axis(pGazeEpoch.scene.limit);
            axis('equal');
                                                           
            % keypress call back     
            fprintf('Please select saccade: ("k",previous, "j":next saccade, "y":mark select, "n":mark unselect,"q":quit, "r":reload)\n');
            set(h.fig,'KeyPressFcn',@keypress_callback)
            set(h.fig,'CloseRequestFcn',@closeRequest_callback);
            waitfor(h.fig)
    end
    
    function keypress_callback(src,evt)
        local_activeInd       = src.UserData.activeInd;
        local_flagSelected    = src.UserData.flagSelected;
        local_flagProcessed   = src.UserData.flagProcessed;
        switch evt.Key
            case 'k' % go to previous not processed saccade
                notProcessed = find(~local_flagProcessed);
                indInNotProcessed = find(notProcessed==local_activeInd);
                if indInNotProcessed == 1
                    newSaccadeInd = notProcessed(1);
                else
                    newSaccadeInd = notProcessed(indInNotProcessed-1);
                end
                src.UserData.activeInd =  newSaccadeInd;
                updatefigure(src);                       
            case 'j' % go to next not processed saccade
                notProcessed = find(~local_flagProcessed);
                indInNotProcessed = find(notProcessed==local_activeInd);
                if indInNotProcessed == numel(notProcessed)
                    newSaccadeInd = notProcessed(end);
                else
                    newSaccadeInd = notProcessed(indInNotProcessed+1);
                end
                src.UserData.activeInd =  newSaccadeInd;
                updatefigure(src);                      
            case 'y' % mark select and go to next unprocessed saccades
                local_flagSelected(local_activeInd) = 1;
                local_flagProcessed(local_activeInd) = 1;
                newSaccadeInd = find(~local_flagProcessed,1);                                              
                src.UserData.activeInd =  newSaccadeInd;
                src.UserData.flagSelected = local_flagSelected;
                src.UserData.flagProcessed = local_flagProcessed;               
                if sum(local_flagSelected) == n
                    fprintf('selection completed, %d saccades selected\n',sum(src.UserData.flagSelected))
                    close(src);return;
                end               
                if isempty(newSaccadeInd)
                    fprintf('all processed, %d saccades selected\n',sum(src.UserData.flagSelected))
                    close(src);return;
                end
                updatefigure(src);                         
            case 'n' % mark unselect and go to next unprocessed saccade
                local_flagProcessed(local_activeInd) = 1;
                newSaccadeInd = find(~local_flagProcessed,1);
                src.UserData.activeInd =  newSaccadeInd;
                src.UserData.flagProcessed = local_flagProcessed;                
                if isempty(newSaccadeInd)
                    fprintf('all processed!\n')
                    close(src);return;
                end
                updatefigure(src);                                                 
            case 'r' % reload
                src.UserData.activeInd = 1;
                src.UserData.flagSelected = zeros(size(src.UserData.flagSelected));
                src.UserData.flagProcessed = zeros(size(src.UserData.flagProcessed));
                updatefigure(src);                                                 
            case 'q'% quit and confirm current selection
                fprintf('advanced quit, %d saccades selected\n',sum(src.UserData.flagSelected));
                close(src);return;
            otherwise
                warning('("k",previous, "j":next saccade, "y":mark select, "n":mark unselect,"q":quit, "r":reload)\n')
        end               
    end
    
    function closeRequest_callback(src,~)
        selected_saccades = src.UserData.saccades(logical(src.UserData.flagSelected));
        delete(gcf);
    end

    function updatefigure(callbackObj)
        foo = callbackObj.UserData;
        set(h.fig,'Name',sprintf('saccade %03d/%03d, %03d processed, %03d selected',...
            foo.activeInd,numel(foo.saccades), sum(foo.flagProcessed), sum(foo.flagSelected)))
        set(h.xtrace,{'XData','yData'},{foo.saccades(foo.activeInd).time, foo.saccades(foo.activeInd).position(:,1)});
        set(h.ytrace,{'XData','yData'},{foo.saccades(foo.activeInd).time, foo.saccades(foo.activeInd).position(:,2)});
        set(h.trajectory,{'xData','yData'},{foo.saccades(foo.activeInd).position(:,1),foo.saccades(foo.activeInd).position(:,2)});
    end
end

