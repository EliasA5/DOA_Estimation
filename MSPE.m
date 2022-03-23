function [err] = MSPE(theta, theta_est)
    err = (cos(theta) - cos(theta_est)).^2 + (sin(theta) - sin(theta_est)).^2;
end