function [angle_output] = wrap_2pi(angle)
 angle = mod(angle + pi, 2*pi);
    if angle < 0
        angle = angle + 2*pi;
    end

    angle_output = angle - pi;
end

