function out = generic_kalman(cfg_in,x)
% function out = generic_kalman(cfg_in,x)
%

cfg_def.sigma_model = 10^5;
cfg_def.P_model = 10^-3;
cfg_def.sigma_meas = 50;
cfg_def.dt = 1/60;

cfg = ProcessConfig(cfg_def,cfg_in);


nP = length(x);

% initial state
Xk_prev = [x(1); x(60)-x(1)]; % assumes 60 samples/s

% current state
Xk = [];

% motion equation
Phi = [1 cfg.dt;
       0  1];

% error/confidence matrix
P = [cfg.P_model 0;
     0 cfg.P_model];

% process noise covariance
Q = [(cfg.dt^4)/4  (cfg.dt^3)/3;
     (cfg.dt^3)/3  cfg.dt^2] * cfg.sigma_model;
             
% measurement matrix -- only measure position
M = [1 0];

% measurement noise covariance
R = cfg.sigma_meas^2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Kalman iteration %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Buffers for later display
Xk_buffer = zeros(2,nP);
Xk_buffer(:,1) = Xk_prev;
Z_buffer = zeros(1,nP);

for k = 1:nP-1
    
    % Z is the measurement vector. In our
    % case, Z = TrueData + RandomGaussianNoise
    Z = x(k+1); %+sigma_meas*randn;
    Z_buffer(k+1) = Z;
    
    % Kalman iteration
    P1 = Phi*P*Phi' + Q;
    S = M*P1*M' + R;
    
    % K is Kalman gain. If K is large, more weight goes to the measurement.
    % If K is low, more weight goes to the model prediction.
    K = P1*M'*inv(S);
    
    % if missing data, set K = 0
    if isnan(Z)
        K = zeros(size(K));
        Z = 0; % to prevent NaNs
    end
    
    P = P1 - K*M*P1;
    
    Xk = Phi*Xk_prev + K*(Z-M*Phi*Xk_prev);
    Xk_buffer(:,k+1) = Xk;
    
    % For the next iteration
    Xk_prev = Xk; 
end;

out = Xk_buffer(1,:);