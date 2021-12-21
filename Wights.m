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
figure;
subplot(1,2,1)
imagesc(log10(abs(W)));
W( W < c) = 0;  

subplot(1,2,2)
imagesc(log10(abs(W)));
save('Weights.mat',"W");
writematrix(W,'Weights.txt');

D = diag(sum(W,2));   % Diagonal matrix that hold the sum of all the weights for each edge
L = D - W;    % Non-normalized Laplacian matrix

