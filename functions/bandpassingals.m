function sig=bandpassingals(signal,fs)

[b,a]=butter(2,[20 500]/(fs/2));
sig=filtfilt(b,a,signal.').';

end