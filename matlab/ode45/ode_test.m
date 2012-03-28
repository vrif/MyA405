function ode_test
  % test prograrm to demonstrate how to
  % use an event handler to stop an
  % integration and savethe results
  % do "help ode45" and scan for EVENTS
  
  tspan=0:0.01:10;
  [t,y]=ode45(@F,tspan,[2,8]);
  figure(1);
  clf;
  plot(t,y(:,1));
  title('height vs. time, no stopping');

  stopHeight=5;
  options=odeset('Events',@events);
  [t,y]=ode45(@F,tspan,[2,8],options);
  figure(2);
  clf;
  plot(t,y(:,1));
  title('height vs. time, stop when descending past 5');

  function [value,isterminal,direction] = events(t,y)
        % Locate the time when height passes through stopHeight in a decreasing direction
        % and stop integration.  Here we use a nested function to avoid
        % passing the additional parameter stopHeight as an input argument.
        value = y(1) - stopHeight;     % detect height = 0
        isterminal = 1;   % stop the integration
        direction = -1;   % only when approaching stopHeight from above
  end

  
end

function xp=F(t,x)
  xp=zeros(2,1); % since output must be a column vector
  xp(1)=x(2);
  xp(2)= -t*x(1) - exp(t)*x(2) + 3*sin(2*t);
end
