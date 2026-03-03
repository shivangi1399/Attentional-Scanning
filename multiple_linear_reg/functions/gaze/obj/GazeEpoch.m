classdef GazeEpoch < handle
    % Eyedata Summarizer class
    properties (SetAccess=protected)
        time       % nx1 double, in seconds
        position   % nx2 double, [x,y], in pix
        velocity   % nx2 double, [x,y], in degree/sec
        pixPerDeg   % pixels per degree 
    end
    
    properties (SetAccess=protected, Hidden)
        scene = struct('display',repmat(0.5*ones(1000),1,1,3),...
            'limit', [-500, 500, -500, 500]);
        parent      % handle obj
        children    % handle array, [sac01,sac02]
    end
       
    methods
        function obj = GazeEpoch(position,time,pixPerDeg)
            if nargin == 0; return; end % to support array init
            obj.position = position;
            obj.time = time;
            obj.pixPerDeg = pixPerDeg;
            obj.velocity = calculaeVelocity;
            
            function velocity = calculaeVelocity
                foo = obj.position(5:end,:)+obj.position(4:end-1,:) ...
                    -obj.position(2:end-3,:)-obj.position(1:end-4,:);
                dt  = mean(diff(obj.time));
                velocity = ([nan(2,2);foo./(6*dt);nan(2,2)])./obj.pixPerDeg;
            end  
            
        end             
        function setPixPerDeg(obj,ppd)
            obj.pixPerDeg = ppd;
        end
        function setParent(obj,parent)
            obj.parent = parent;
        end
        function setChildren(obj,children)
            obj.children = children;
        end        
        function addChild(obj,child)
            obj.children(end+1) = child;
        end
        % scene methods
        function showScene(obj)
            image(obj.scene.limit(1:2),obj.scene.limit(4:-1:3),...
                    obj.scene.display);
            set(gca,'YDir','normal');
            set(gca,'Visible','off');
            axis(obj.scene.limit);
            axis('equal');
        end
        function setSceneDisplay(obj,display)
            obj.scene.display = display;
        end
        function setSceneLimit(obj,limit)
            assert(isa(limit,'double'))
            assert(numel(limit)==4)
            obj.scene.limit = limit;
        end
        function clearScene(obj)
            obj.scene.disaply = repmat(0.5*ones(1000),1,1,3);
            obj.scene.limit = [-500,500,-500,500];
        end
        
    end        
end