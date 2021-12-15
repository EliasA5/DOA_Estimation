clear all
close all
clc

%------------------------------------------------------
data = load("L_A.mat");
L_A = data.weights;
latitude = L_A(:,1);
longtitude = L_A(:,2);
R = 6371;   % radius of the earth = 6371 km
X = R*cosd(latitude).*cosd(longtitude);
Y = R*cosd(latitude).*sind(longtitude);
Z = R*sind(latitude);
W = zeros(length(X));
d = length(X);    % check what those parameters mean or should be
sigma = 1;
c = 10^(-5);    % the threshold is 10^(-7)

for i=1:length(X)
    W(i,:) = exp(-d*sqrt((X(i)-X).^2 + (Y(i)-Y).^2 + (Z(i)-Z).^2)/(2*sigma));
end

W( W < c) = 0;  
save('Weights.mat',"W");

L = diag(W) - W;

