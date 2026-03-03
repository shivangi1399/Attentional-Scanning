function plotSaccades(saccades, tlim)
%PLOT_SACCADE Show saccade positon and velocity profile    
% 
%   PLOT_SACCADE(saccades) plot saccades with parent GazeEpoch extend 0.1s
%   before(after) saccade start(end).
%
%   PLOT_SACCADE(saccades) plot saccades with parent GazeEpoch from
%   saccadeStart + tlim(1) to saccadeEnd + tlim(2)
%   Example:
%   plot_saccade(saccades, [-0.1,0.1]

assert(isa(saccades,'Saccade'), 'Not a Saccade Instance')

if nargin == 1
    tlim = [-0.1,0.1];
elseif nargin == 2
    assert(isa(tlim,double),'Invalid datatype for tlim');
    assert(numel(tlim)==2, 'Invalid numel for tlim');
end

n = numel(saccades);
h = cell(n,1);

for ii = 1 : numel(saccades)
    thisSaccade = saccades(ii);
    thisParentGaze = saccades(ii).parent;
    startIndInParent = time2ind(thisSaccade.time(1)+tlim(1),thisParentGaze.time);
    endIndInParent = time2ind(thisSaccade.time(end)+tlim(2),thisParentGaze.time);
    
    h{ii}.fig = figure(ii);
    set(h{ii}.fig,'Name', sprintf('Saccade Number %d out of %d',ii,n));
    % xtrace
    subplot(2,4,1)
    xdata = thisParentGaze.time(startIndInParent:endIndInParent);
    ydata = thisParentGaze.position(startIndInParent:endIndInParent,1)./thisParentGaze.pixPerDeg;    
    plot(xdata,ydata,'k-'); hold on;
    
    xdata = thisSaccade.time;
    ydata = thisSaccade.position(:,1)./thisSaccade.pixPerDeg;    
    h{ii}.xtrace = plot(xdata,ydata,'r-');
    title('Horizontal component')
    ylabel(['Position [',char(176),']']);

    % xvelocity
    subplot(2,4,5)  
    xdata = thisParentGaze.time(startIndInParent:endIndInParent);
    ydata = thisParentGaze.velocity(startIndInParent:endIndInParent,1);

    plot(xdata,ydata,'k-'); hold on;
    
    xdata = thisSaccade.time;
    ydata = thisSaccade.velocity(:,1);
    h{ii}.xvel = plot(xdata,ydata,'r-');
    xlabel('Time [sec]');
    ylabel(['Velocity [',char(176),']']);
    
    % ytrace
    subplot(2,4,2)
    xdata = thisParentGaze.time(startIndInParent:endIndInParent);
    ydata = thisParentGaze.position(startIndInParent:endIndInParent,2)./thisParentGaze.pixPerDeg;   
    plot(xdata,ydata,'k-'); hold on;
    
    xdata = thisSaccade.time;
    ydata = thisSaccade.position(:,2)./thisSaccade.pixPerDeg;
    h{ii}.ytrace = plot(xdata,ydata,'r-');
    title('Vertical component')
    
    % yvelocity
    subplot(2,4,6)
    xdata = thisParentGaze.time(startIndInParent:endIndInParent);
    ydata = thisParentGaze.velocity(startIndInParent:endIndInParent,2);   
    plot(xdata,ydata,'k-'); hold on;
    
    xdata = thisSaccade.time;
    ydata = thisSaccade.velocity(:,2);       
    h{ii}.yvel = plot(xdata,ydata,'r-');
    xlabel('Time [sec]');
    
    % 2d trajectory
    subplot(2,4,[3 4 7 8])
    image(thisParentGaze.scene.limit(1:2),thisParentGaze.scene.limit(4:-1:3),...
            thisParentGaze.scene.display);hold on;
    set(gca,'YDir','normal');
    set(gca,'Visible','off');
    plot(thisParentGaze.position(startIndInParent:endIndInParent,1),thisParentGaze.position(startIndInParent:endIndInParent,2),'-','Color',[0.25,0.25,0.25],'LineWidth',0.5);
    h{ii}.trajectory = plot(thisSaccade.position(:,1),thisSaccade.position(:,2),'.-r','LineWidth',1);
    axis(thisParentGaze.scene.limit);
    axis('equal');
 end