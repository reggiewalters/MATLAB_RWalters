% This function performs an ensemble reordering method called the Schaake
% Shuffle (Clark et al., 2004) on an ensemble of forecasts (X), given an
% identically-sized ensemble of corresponding historical observations (Y)
% for the same or very similar day(s) of year.
% Input dimensions are [M x N], where M = number of ensemble members and 
% N = forecast days.
% R. Walters, HHWP, Nov 2019

function X_ss = Schaake_Shuffle(X, Y)

% % % Sort forecast input array X
Chi = sort(X);

% % % Sort observed input aray Y and get index vector B
[~, B] = sort(Y);

% % % Perform Schaake Shuffle for the re-ordered array of forecast ensemble
%     members
[n, L] = size(X);           % number of ensemble members, forecast length

X_ss = nan.*ones(n,L);
for i = 1:L
    [~, ~, q] = intersect(repmat((1:n)', 1, n), B(:,i));
    X_ss(:,i) = Chi(q,i);
end
