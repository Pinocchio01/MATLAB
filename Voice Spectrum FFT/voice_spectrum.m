[x1, fs] = audioread('Lost Cause.mp3'); % data and sampling frequency
x1 = 0.5*(x1(:,1) + x1(:,2)); % average of two tracks

subplot(3,1,1);
plot(x1);
title('Original voice signal');
xlabel('Time');
ylabel('Volume');

y1 = fft(x1);   % (0, 2*pi)
% plot(abs(y1)); % original spectrum, y1 complex number
y1 = fftshift(y1); % cut off the spectrum in (pi, 2*pi) and move to front (-pi, 0), basically a time shift property of Fourier Transform


% FFT calculate vertical values but no horizontal values, below calculate

derta_fs = fs/length(x1);

% real function, amplitude specturm even function, phase spectrum odd
% function
subplot(3,1,2);
% pay attention to discrete value FT and continuous FT (horizontal *T, vertical /T)
plot(-fs/2:derta_fs:fs/2-derta_fs, abs(y1)/fs);   % in Hz, max w_M = fs/2 = 22050 Hz
title('Original voice frequency amplitude spectrum');


subplot(3,1,3);
plot(-fs/2:derta_fs:fs/2-derta_fs, atan2(imag(y1), real(y1)));
title('Original voice frequency phase spectrum');