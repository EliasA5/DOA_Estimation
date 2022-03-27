function [err] = MSPE(theta, theta_est, type)
switch(type)
    case 'cyclic'
        err = (cos(theta) - cos(theta_est)).^2 + (sin(theta) - sin(theta_est)).^2;
    case 'MSPE'
        err = custom_mod2pi(theta-theta_est).^2;
    case 'MSE'
        err = (theta - theta_est).^2;
end