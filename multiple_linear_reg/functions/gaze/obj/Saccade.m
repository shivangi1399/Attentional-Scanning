classdef Saccade < GazeEpoch        
    properties(Dependent)        
        magnitude   % 1x1 double, in degree
        direction   % 1x1 double, in radians
        duration    % 1x1 double, in seconds
        peakVel     % 1x1 double, in degree/sec
    end
    
    properties(Dependent, Hidden)
        vector      % 1x1 double, complex
    end
    
    properties(SetAccess=private, Hidden)
        info
        startIndInParent
        endIndInParent
    end
    
    methods
        function obj = Saccade(position,time,pixPerDeg)
            if nargin == 0
                superArgs = {};
            elseif nargin == 2
                superArgs = {position, time, pixPerDeg};
            else
                error('Wrong number of input arguments')
            end
            obj@GazeEpoch(superArgs{:});
        end        
        function fetchDataFromParent(obj,startInd,endInd)
            obj.startIndInParent = startInd;
            obj.endIndInParent = endInd;
            obj.position = obj.parent.position(startInd:endInd,:);
            obj.time = obj.parent.time(startInd:endInd,:);
            obj.velocity = obj.parent.velocity(startInd:endInd,:);
            obj.scene = obj.parent.scene;
            obj.pixPerDeg = obj.parent.pixPerDeg;
        end 
        function logInfo(obj,info)
            obj.info = info;
        end
        
        % get methods
        function vector = get.vector(obj)
            assert(~isempty(obj.position),'Saccade Object property `position` is empty')
            foo = obj.position(end,:) - obj.position(1,:);
            vector = (foo(1) + 1i*foo(2))./obj.pixPerDeg;
        end
        function magnitude = get.magnitude(obj)
            magnitude = abs(obj.vector);
        end        
        function direction = get.direction(obj)
            direction = angle(obj.vector);
        end
        function duration = get.duration(obj)
            assert(~isempty(obj.time),'Saccade Object property `time` is empty');
            duration = obj.time(end)-obj.time(1);            
        end
        function peakVel = get.peakVel(obj)
            assert(~isempty(obj.velocity),'Saccade Object property `velocity` is empty')
            peakVel = max(abs(obj.velocity));
        end
        
        function trace = normTrace(obj)
            complexPos = complex(obj.position(:,1), obj.position(:,2));
            complexPosNorm = (complexPos-complexPos(1)) ./ (complexPos(end)-complexPos(1)); 
            x = real(complexPosNorm);
            y = imag(complexPosNorm);
            trace = [x,y];
        end
        
        % visualization methods, raw data
        function plot(obj, varargin)
        % PLOT(obj,options) plot saccade (with options):
        %   options: cell array of strings specifying:
        %       'noParent': don't show parent gazeEpoch
        %       'noScene': don't show associated scene      
            allowedOptions = {
                'noParent',...
                'noScene',...
                'fullScene'
            };
        
            if nargin == 1; varargin = {};end
            if ~isempty(varargin)
                assert(all(ismember(varargin,allowedOptions)),...
                    'Not valid option');
            end
            
            flag_noParent = any(cellfun(@(x) strcmp(x,'noParent'),varargin));
            flag_noScene = any(cellfun(@(x) strcmp(x,'noScene'),varargin));
            flag_fullScene = any(cellfun(@(x) strcmp(x,'fullScene'),varargin));
        
            thisSaccade = obj;
            thisParentGaze = obj.parent;
            
            % plot scene
            set(gcf,'Name','Saccade Trajectory');
            scene4plot = thisParentGaze.scene;

            if flag_noScene
                sceneAlpha = 0;
            else
                sceneAlpha = 1;
            end
            
            if flag_fullScene
                image(scene4plot.limit(1:2),scene4plot.limit(4:-1:3),...
                    scene4plot.display,...
                    'AlphaData',sceneAlpha);
            else % center the saccade
                saccadeStartXY  = thisSaccade.position(1,:);
                viewSizeInPix   = thisSaccade.magnitude * thisSaccade.pixPerDeg+50;
                image([saccadeStartXY(1)-viewSizeInPix,saccadeStartXY(1)+viewSizeInPix],...
                      [saccadeStartXY(2)+viewSizeInPix,saccadeStartXY(2)-viewSizeInPix],...
                    scene4plot.display,...
                    'AlphaData',sceneAlpha);
            end
            
            hold on;
            set(gca,'YDir','normal');
            set(gca,'Visible','off');
            
            % plot parent gaze epoch
            if ~flag_noParent
                if ~isempty(obj.parent)
                    plot(thisParentGaze.position(:,1),...
                        thisParentGaze.position(:,2),...
                        '-','Color',[0.7,0.7,0.7],'LineWidth',0.5);
                end
            end
            
            % plot current saccade
            xdata = thisSaccade.position(:,1)';
            ydata = thisSaccade.position(:,2)';
            zdata = zeros(size(xdata));
            cdata = 1:numel(xdata);
            surf([xdata;xdata],[ydata;ydata],[zdata;zdata],[cdata;cdata],...
                'FaceColor','none',...
                'EdgeColor','interp',...
                'LineWidth',1);            
            axis('equal');
            scatter(xdata,ydata,10,cdata,'filled');
            colormap(flip(autumn));
            hold off;
        end        
        function plot_positionProfile(obj, tBefore, tAfter)
            assert(nargin == 1 | nargin == 3, 'Wrong number of input'); 
            if nargin == 1 % use default time offset
                tBefore = 0.010;
                tAfter = 0.010;
            end
                          
            thisSaccade = obj;
            try
                assert(~isempty(thisSaccade.parent))
                thisParentGaze = obj.parent;
            catch
                % no parent gazeEpoch available, use itself as parent
                thisParentGaze = thisSaccade;            
            end
            
            ind1 = time2ind(thisSaccade.time(1)-tBefore, thisParentGaze.time);  %plot start
            ind2 = thisSaccade.startIndInParent;    % saccade start
            ind3 = thisSaccade.endIndInParent;      % saccade end
            ind4 = time2ind(thisSaccade.time(end)+tAfter,thisParentGaze.time); %plot end
            
            tdata = transpose((thisParentGaze.time(ind1:ind4)-thisSaccade.time(1))*1000); % in ms, rel. sac onset
            xdata = transpose(thisParentGaze.position(ind1:ind4,1)/obj.pixPerDeg); % in degree
            ydata = transpose(thisParentGaze.position(ind1:ind4,2)/obj.pixPerDeg); % in degree
            zdata = zeros(size(tdata));
            cdata = [ones(1,ind2-ind1+1),1:(ind3-ind2),ones(1,ind4-ind3)];
            
            % x position profile
            set(gcf,'Name','Saccade Position Profile')
            subplot(2,1,1)
            scatter(tdata,xdata,10,cdata,'filled');
            hold on;
            surf([tdata;tdata],[xdata;xdata],[zdata;zdata],[cdata;cdata],...
                'FaceColor','none',...
                'EdgeColor','interp',...
                'LineWidth',1);
            colormap(flip(autumn));
            plot([0,0],ylim,'--k');
            plot([thisSaccade.duration*1000,thisSaccade.duration*1000],ylim,'--k');
            title('X-Position');
            ylabel('\DeltaPosition [dva]')
            
            % y position profile
            subplot(2,1,2)
            scatter(tdata,ydata,10,cdata,'filled');
            hold on;
            surf([tdata;tdata],[ydata;ydata],[zdata;zdata],[cdata;cdata],...
                'FaceColor','none',...
                'EdgeColor','interp',...
                'LineWidth',1);
            colormap(flip(autumn));
            plot([0,0],ylim,'--k');
            plot([thisSaccade.duration*1000,thisSaccade.duration*1000],ylim,'--k');
            title('Y-Position');
            xlabel('Time relative to saccade onset [ms]');
            ylabel('\DeltaPosition [dva]')
        end
        function plot_velocityProfile(obj, tBefore, tAfter)
            assert(nargin == 1 | nargin == 3,'Wrong number of input');
            if nargin == 1 % use default time offset
                tBefore = 0.010;
                tAfter = 0.010;
            end
            
            thisSaccade = obj;
            try
                assert(~isempty(thisSaccade.parent))
                thisParentGaze = obj.parent;
            catch
                % no parent gazeEpoch available, use itself as parent
                thisParentGaze = thisSaccade;            
            end
            
            ind1 = time2ind(thisSaccade.time(1)-tBefore, thisParentGaze.time);  %plot start
            ind2 = thisSaccade.startIndInParent;    % saccade start
            ind3 = thisSaccade.endIndInParent;      % saccade end
            ind4 = time2ind(thisSaccade.time(end)+tAfter,thisParentGaze.time); %plot end
            
            tdata = transpose((thisParentGaze.time(ind1:ind4)-thisSaccade.time(1))*1000); % in ms, rel. sac onset
            xdata = transpose(thisParentGaze.velocity(ind1:ind4,1)); % in degree/sec
            ydata = transpose(thisParentGaze.velocity(ind1:ind4,2)); % in degree/sec
            zdata = zeros(size(tdata));
            cdata = [ones(1,ind2-ind1+1),1:(ind3-ind2),ones(1,ind4-ind3)];
            
            % x velocity profile
            set(gcf,'Name','Saccade Velocity Profile')
            subplot(2,1,1)
            scatter(tdata,xdata,10,cdata,'filled');
            hold on;
            surf([tdata;tdata],[xdata;xdata],[zdata;zdata],[cdata;cdata],...
                'FaceColor','none',...
                'EdgeColor','interp',...
                'LineWidth',1);
            colormap(flip(autumn));
            plot([0,0],ylim,'--k');
            plot([thisSaccade.duration*1000,thisSaccade.duration*1000],ylim,'--k');
            title('X-Velocity');
            ylabel('Velocity [dva/sec]')
            
            % y velocity profile
            subplot(2,1,2)
            scatter(tdata,ydata,10,cdata,'filled');
            hold on;
            surf([tdata;tdata],[ydata;ydata],[zdata;zdata],[cdata;cdata],...
                'FaceColor','none',...
                'EdgeColor','interp',...
                'LineWidth',1);
            colormap(flip(autumn));
            plot([0,0],ylim,'--k');
            plot([thisSaccade.duration*1000,thisSaccade.duration*1000],ylim,'--k');
            title('Y-Velocity');
            xlabel('Time relative to saccade onset [ms]');
            ylabel('Velocity [dva/sec]')
            
        end       
        function plot_2dVelocityProfile(obj, tBefore, tAfter)
            assert(nargin == 1 | nargin == 3,'Wrong number of input');
            if nargin == 1 % use default time offset
                tBefore = 0.010;
                tAfter = 0.010;
            end
            
            thisSaccade = obj;
            try
                assert(~isempty(thisSaccade.parent))
                thisParentGaze = obj.parent;
            catch
                % no parent gazeEpoch available, use itself as parent
                thisParentGaze = thisSaccade;            
            end
            
            %plot decision boundary (elliptic)
            fcn_plot_ellipse = @(a, b, lineWidth) ...
                plot(a*cos(linspace(0,2*pi,100)),b*sin(linspace(0,2*pi,100)),'--k',...
                    'LineWidth',lineWidth);
            
            v = thisParentGaze.velocity;
            v_median = median(v,1,'omitnan');
            v_sqr = v .^ 2;
            v_sqr_median = median(v_sqr,1,'omitnan');
            v_sigma = sqrt(v_sqr_median - v_median.^2);
            
            set(gcf,'Name','2D-Velocity Profile');
            eta = 3*v_sigma;
            fcn_plot_ellipse(eta(1),eta(2),1.5); hold on;
            eta = 6*v_sigma;
            fcn_plot_ellipse(eta(1),eta(2),1.0);
            eta = 10*v_sigma;
            fcn_plot_ellipse(eta(1),eta(2),0.5);
            legend({'\lambda=3','\lambda=6','\lambda=10'},...
                'Location', 'SouthWest');
            
            %plot velocity traces
            ind1 = time2ind(thisSaccade.time(1)-tBefore, thisParentGaze.time);  %plot start
            ind2 = thisSaccade.startIndInParent;    % saccade start
            ind3 = thisSaccade.endIndInParent;      % saccade end
            ind4 = time2ind(thisSaccade.time(end)+tAfter,thisParentGaze.time); %plot end
            
            xdata = transpose(thisParentGaze.velocity(ind1:ind4,1)); % in degree/sec
            ydata = transpose(thisParentGaze.velocity(ind1:ind4,2)); % in degree/sec
            zdata = zeros(size(xdata));
            cdata = [ones(1,ind2-ind1+1),1:(ind3-ind2),ones(1,ind4-ind3)];
            
            scatter(xdata,ydata,10,cdata,'filled','HandleVisibility','off');
            surf([xdata;xdata],[ydata;ydata],[zdata;zdata],[cdata;cdata],...
                'FaceColor','none',...
                'EdgeColor','interp',...
                'LineWidth',1,...
                'HandleVisibility','off');
            colormap(flip(autumn));
            title('2D-Velocity');
            axis('equal');
            xlabel('X-Velocity [dva/sec]');
            ylabel('Y-Velocity [dva/sec]');
        end
        
        % visualization methods, normalized
        function plot_norm(obj)
        %Plot normalized trajectory
           complexPos = complex(obj.position(:,1), obj.position(:,2));
           complexPosNorm = (complexPos-complexPos(1)) ./ (complexPos(end)-complexPos(1));
           
           xdata = real(complexPosNorm)';
           ydata = imag(complexPosNorm)';
           zdata = zeros(size(xdata));
           cdata = 1:numel(xdata);
           
           scatter(xdata,ydata,10,cdata,'filled'); hold on;
           surf([xdata;xdata],[ydata;ydata],[zdata;zdata],[cdata;cdata],...
            'FaceColor','none',...
            'EdgeColor','interp',...
            'LineWidth',1);
           axis('equal');
           colormap(flip(autumn));
        end
    end
end