function time = chutetime(alt_range, ref_rate, opts)
    % Calculate descent time through altitude range
    % time = chutetime(alt_range, ref_rate, [, Name, Value])
    %   Inputs
    %   alt_range   (double)    Range of altitudes to descend through
    %                           [start_alt] or [start_alt end_alt], end_alt defaults to 0
    %   ref_rate    (double)    Reference descent rate 
    %   Name-value pairs
    %   ref         ({0} | double)  Altitude at which ref_rate is calculated (default 0)
    %   ground      ({0} | double)  Ground-level altitude added to *all* altitude inputs
    %   units       ({"m"} | "ft")  Length units for altitudes and rate
    arguments
        alt_range (1,:) double {mustBeNonnegative};
        ref_rate (1,1) double {mustBePositive};
        opts.ref (1,1) double {mustBeNonnegative} = 0;
        opts.units (1,1) string {mustBeMember(opts.units, ["m", "ft"])} = "m";
        opts.ground (1,1) double {mustBeNonnegative} = 0;
    end

    % populate input defaults
    if isscalar(alt_range)
        alt_range = [alt_range 0];
    end
    opts.ref = opts.ref + opts.ground;
    alt_range = alt_range + opts.ground;

    m_per_ft = 1/3.048;

    switch opts.units
        case "m"
        case "ft"
            ref_rate = ref_rate * m_per_ft;
            opts.ref = opts.ref * m_per_ft;
            alt_range = alt_range * m_per_ft;
        otherwise 
            error("Unrecognized unit");
    end

    start_height = max(alt_range);
    end_height = min(alt_range);
    ref_density = density_aloft(opts.ref);

    % calculate maximum possible time using slowest possible descent rate
    % in reality, ode45() will terminate using @finished_descent
    min_rate = ref_rate * (ref_density / density_aloft(end_height));
    max_time = (start_height - end_height) / min_rate;

    % use ode45 to numerically find final time
    options = odeset(Events = @finished_descent, AbsTol = 1e-10, RelTol = 1e-10);
    [time, hgt] = ode45(@descent_rate, [0 max_time], start_height, options);

    % check that output is sane (within 10m of actual)
    assert(abs(hgt(end) - end_height) <= 10);

    time = time(end);

    % descent ODE
    function rate = descent_rate(~, h)
        rate = -ref_rate * sqrt(ref_density / density_aloft(h));
    end

    % event function indicating finished descent
    function [zeroing, isterminal, direction] = finished_descent(~, h)
        zeroing = h - end_height;
        isterminal = 1;
        direction = -1;
    end
end

function rho = density_aloft(alt_asl_m)
    if (alt_asl_m <= 11e3)
        rho = 1.225 * (1 - 22.558e-6 * alt_asl_m) ^ 4.2559;
    elseif (11e3 < alt_asl_m) && (alt_asl_m <= 20e3)
        rho = 0.3639 * exp(-157.69e-6 * (alt_asl_m - 11000));
    else
        error("Input out of range");
    end
end

