function [val] = steeringVector(rm, theta, alpha, lambda) %as described in https://www.antenna-theory.com/definitions/steering.php
    val = exp(-1i .* 2*pi/lambda*rm*[sin(theta)*cos(alpha), sin(theta)*sin(alpha), cos(theta)].');
end