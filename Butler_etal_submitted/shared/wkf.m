function [z,C] = wkf(y,A,Cx,cy,x0,C0)
%% Wrapped Kalman filter
% function [z,C] = wkf(y,A,Cx,cy,x0,C0)
%
% Inputs: y  - Tx1 obs sequence
%         A  - NxN state transition matrix
%         Cx - NxN state covariance
%         cy - obs variance
%         x0 - Nx1 initial state (default: zeros(N,1))
%         C0 - NxN initial state covariance (default: Cx)
% Outputs: z - NxT state sequence estimates
%          C - NxN estimated state covariance (converges)
%
% The wrapped Kalman filter infers the hidden state sequence of a LDS
% where the state lies on the unit circle. The state position is modeled
% with a wrapped Gaussian (wG) and all other components (2 to N) are
% modeled with Gaussians. The observation is a noisy measurement of
% the state position.
%
% Algorithm is described in:
%    "A Wrapped Kalman Filter for Azimuthal Speaker Tracking"
%    IEEE Signal Processing Letters, December 2013, Volume: 20 Issue: 12
%    Pages 1257-1260 (Johannes Traa and Paris Smaragdis)
%
% Johannes Traa - UIUC 2013

T = length(y);
N = size(A,1);

%% check inputs
if nargin < 5 || isempty(x0); x0 = zeros(N,1); end
if nargin < 6 || isempty(C0); C0 = Cx;         end

%% params
B = [1 zeros(1,N-1)]; % observation matrix
L = 1; % only 3 terms necessary
ll = 2*pi*(-L:L)'; % wG terms

%% initialize
z = zeros(N,T);
z(:,1) = x0;
C = C0;

%% filter
for t=2:T
  % predict
  z(:,t) = A*z(:,t-1);
  C = A*C*A' + Cx;
  
  % posteriors
  %yp = B*z(:,t) + ll; % predicted obs
  yp = B*z(:,t);
  %yt = round(yp(L+1)/(2*pi))*2*pi + y(t); % compensate for wG drift
  yt = y(t);
  %pp = -(diffang_rad(yt,yp).^2)/(2*cy);
  %pp = exp(bsxfun(@minus,pp,logsum(pp,1)));
  
  % correct
  K = (C*B')/(B*C*B' + cy);
  %yy = pp'*(diffang_rad(yt,yp)); % weighted combination of innovations
  yy = diffang_rad(yt,yp); % weighted combination of innovations
  
  % if missing data, just go with prediction
  if isnan(yt)
      z(:,t) = A*z(:,t-1);
      z(1,t) = add2pi(z(1,t),0);
  else
      kalman_update = K*yy;
      z(1,t) = add2pi(z(1,t),kalman_update(1));
      z(2:end,t) = kalman_update(2:end);
  end
    

  
  C = (eye(N) - K*B)*C;
end







function y = logsum(x,d)
%% Log-sum-exp for sumation in the log domain
% function y = logsum(x,d)
%
% Inputs: x - matrix or vector
%         d - dimension to sum over (no default)
% Output: y - logsum output
%
% Johannes Traa - UIUC 2012

m = max(x,[],d);
y = log(sum(exp(bsxfun(@minus,x,m)),d)) + m;