function [value,isterminal,direction] = stopevent(tt,xx,Ap,Aorifice,Cd,rho,speed,S)

if toc > 5 % time in seconds
    value = 0;
    isterminal = 1;     % stop the integration
    direction = 0;      % all events
else
    value = 1;
    isterminal = 0;     % don't stop
    direction = 0;      % all events
end