%%
f1 = 20;
f2 = 25;
f3 = 33.3;

dt = 1e-5; 
t_lag = -0.005; % in s
tvec = -0.5/f1:dt:0.5/f1;

s1 = sin(2*pi*f1*tvec + pi/2);
s1_shifted = sin(2*pi*f1*tvec + 2*pi*t_lag*f1 + pi/2);

s2 = sin(2*pi*f2*tvec + pi/2);
s2_shifted = sin(2*pi*f2*tvec + 2*pi*t_lag*f2 + pi/2);

s3 = sin(2*pi*f3*tvec + pi/2);
s3_shifted = sin(2*pi*f3*tvec + 2*pi*t_lag*f3 + pi/2);

%% 
subplot(321)
plot(tvec*1000,s1,'r');
hold on;
plot(tvec*1000,s1_shifted,'b');
%plot([(1/f1)/4 (1/f1)/4],[-1 1],'k:');
plot([0 0],[-1 1],'k:');
plot([-t_lag*1000 -t_lag*1000],[-1 1],'k:');
set(gca,'YTick',[]); box off; axis tight;

subplot(323)
plot(tvec*1000,s2,'r');
hold on;
plot(tvec*1000,s2_shifted,'b');
%plot([(1/f2)/4 (1/f2)/4],[-1 1],'k:');
plot([0 0],[-1 1],'k:');
plot([-t_lag*1000 -t_lag*1000],[-1 1],'k:');
set(gca,'YTick',[]); box off; axis tight;

subplot(325)
plot(tvec*1000,s3,'r');
hold on;
plot(tvec*1000,s3_shifted,'b');
%plot([(1/f3)/4 (1/f3)/4],[-1 1],'k:');
plot([0 0],[-1 1],'k:');
plot([-t_lag*1000 -t_lag*1000],[-1 1],'k:');
set(gca,'YTick',[]); box off; axis tight;

%%
subplot(322)
fvec = [ 20 25 33.3 ];
phivec = [36 45 60];
[h,ax1,ax2] = plotyy(fvec,phivec,fvec,-phivec);
set(ax1,'Marker','.','MarkerSize',20);
set(ax2,'Marker','.','MarkerSize',20);
hold on;

set(h(1),'YTick',phivec,'XTick',fvec,'YLim',[31 65],'XLim',[18 35.3]);
set(h(2),'YTick',-phivec(end:-1:1),'XTick',fvec,'YLim',[-65 -31],'XLim',[18 35.3]);

[h] = plotyy(fvec,phivec,fvec,-phivec);
hold on;

set(h(1),'YTick',phivec,'XTick',fvec,'YLim',[31 65],'XLim',[18 35.3]);
set(h(2),'YTick',-phivec(end:-1:1),'XTick',fvec,'YLim',[-65 -31],'XLim',[18 35.3]);

grid on;
