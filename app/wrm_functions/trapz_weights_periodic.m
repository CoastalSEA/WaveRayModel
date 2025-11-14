function w = trapz_weights_periodic(theta)
%TRAPZWEIGHTSPERIODIC  Integration weights for periodic trapezoidal rule
%
%   w = trapzWeightsPeriodic(theta)
%
%   theta : N-by-1 vector of angles in radians, assumed increasing and
%           covering [0,2*pi) (last point < 2*pi).
%   w     : N-by-1 vector of trapezoidal weights such that
%           integral ≈ sum(D .* w)
%
%   The weights include the wrap-around interval from theta(end) back to
%   theta(1)+2*pi, so the integration is over the full circle.

    theta = theta(:);
    N = numel(theta);

    % intervals between successive points
    dtheta = diff(theta);
    % closing interval (wrap-around)
    wrap = 2*pi - (theta(end) - theta(1));
    dtheta = [dtheta; wrap];

    % trapezoidal weights: half of left + half of right interval
    w = zeros(N,1);
    for i = 1:N
        left  = dtheta(mod(i-2,N)+1);  % interval before i
        right = dtheta(i);             % interval after i
        w(i)  = (left + right)/2;
    end
end
