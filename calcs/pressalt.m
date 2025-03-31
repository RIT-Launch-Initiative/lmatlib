% Convenience function implementing NOAA's station pressure altitude calculation in various units
% hgt = pressalt(uhgt, press, upress)
%   uhgt:       Output height units, one of "m", "km", "ft", "mi"
%   press:      numeric (array)
%   upress:     Input pressure units, one of "Pa", "kPa", "mbar", "psi"
function hgt = pressalt(uhgt, press, upress)
    arguments
        uhgt (1,1) string {mustBeMember(uhgt, ["m", "ft", "mi", "km"])};
        press double;
        upress (1,1) string {mustBeMember(upress, ["Pa", "kPa", "mbar", "psi"])} = "Pa";
    end
    
    % Convert input
    switch upress
        case "Pa"
            P_mbar = press ./ 100;
        case "kPa"
            P_mbar = press .* 10;
        case "mbar"
            P_mbar = press;
        case "psi"
            P_mbar = 68.94757293168361 .* press;
    end

    hgt_m = (1 - (P_mbar ./ 1013.25).^(0.190284)) .* 145366.45 * 0.3048;

    % Convert output
    switch uhgt
        case "m"
            hgt = hgt_m;
        case "km"
            hgt = hgt_m ./ 1000;
        case "ft"
            hgt = hgt_m .* 3.280839895013123;
        case "mi"
            hgt = hgt_m .* 6.213711922373339e-04;
    end
end

