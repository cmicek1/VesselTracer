function C = ellipsoidInit(epsilon, radius)

sigma = @(u, v, epsilon) ((abs(cos(u) .* cos(v)/radius).^2 + ...
    abs(cos(u) .* sin(v)/radius).^2)...
    .^(1/epsilon) + abs(sin(u)/(radius / 2)).^(2/epsilon)).^(-epsilon/2);

Cx = @(u, v) sigma(u, v, epsilon) .* cos(u) .* cos(v);
Cy = @(u, v) sigma(u, v, epsilon) .* cos(u) .* sin(v);
Cz = @(u, v) sigma(u, v, epsilon) .* sin(u);

set(0,'DefaultFigureVisible','off');
C = fsurf(Cx, Cy, Cz, [-pi/2, pi/2, -pi, pi]);
set(0,'DefaultFigureVisible','on');
end