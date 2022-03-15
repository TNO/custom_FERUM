% Differences between the central points between the elements of x
% at the ends symmetrically extended the last interval
%
%SYNOPSYS
% dc = DIFFC(x)
%
function dc = diffc(x)

if size(x,1) > 1
    x = x';
    type = 'column';
else
    type = 'row';
end

%.....................
% actual calculation
dx  = diff(x);

c   = [x(1)-dx(1)/2, x + [dx/2, dx(end)/2]];
dc  = diff(c);
%.....................

switch type
    case 'column'
        dc = dc';
    case 'row'
        % do nothing
end

end

