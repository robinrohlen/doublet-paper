function [t,F,dFdTc,dFdThr,dFdFmax]=RaikovaForceTwitch5p(fs,Ti,Tlead,Tc,Thr,Ttot,Fmax)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input:  
%  fs      Sample rate (Hz)
%  Ti      Time of stimulus (ms)
%  Tlead   Time delay between stimulus and start of twitch (ms)
%  Thc     Half-contraction time w.r.t. displacement (ms)
%  Tc      Time of contraction, i.e., zero-crossing (ms)
%  Thr     Half-relaxation time w.r.t. displacement (ms)
%  Ttot    Duration of twitch (ms)
%  Fmax    Amplitude of twitch w.r.t. displacement (m/s)
%
% Output:
%  t       Time vector of twitch (ms)
%  F       Tissue velocity vector of twitch (m/s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t=linspace(0,Ttot,Ttot*fs*10^(-3));

F=Fmax.*Tc.^(-log(2.0)./(log(Tc)-log(Thr)+Thr./Tc-1.0)).*t.^(log(2.0)./(log(Tc)-log(Thr)+Thr./Tc-1.0)).*exp(log(2.0)./(log(Tc)-log(Thr)+Thr./Tc-1.0)).*exp(-(t.*log(2.0))./(Tc.*(log(Tc)-log(Thr)+Thr./Tc-1.0)));
dFdTc=Fmax.*exp((Tc.*(log(Tc)-1.0).*log(2.0))./(Tc-Thr+Tc.*log(Thr./Tc))).*exp(((Ti+Tlead-t).*(-log(2.0)))./(Tc-Thr+Tc.*log(Thr./Tc))).*(log(2.0)./(Tc-Thr+Tc.*log(Thr./Tc))+((log(Tc)-1.0).*log(2.0))./(Tc-Thr+Tc.*log(Thr./Tc))-Tc.*log(Thr./Tc).*(log(Tc)-1.0).*1.0./(Tc-Thr+Tc.*log(Thr./Tc)).^2.*log(2.0)).*(-Ti-Tlead+t).^((Tc.*(-log(2.0)))./(Tc-Thr+Tc.*log(Thr./Tc)))-Fmax.*exp((Tc.*(log(Tc)-1.0).*log(2.0))./(Tc-Thr+Tc.*log(Thr./Tc))).*exp(((Ti+Tlead-t).*(-log(2.0)))./(Tc-Thr+Tc.*log(Thr./Tc))).*log(-Ti-Tlead+t).*(log(2.0)./(Tc-Thr+Tc.*log(Thr./Tc))-Tc.*log(Thr./Tc).*1.0./(Tc-Thr+Tc.*log(Thr./Tc)).^2.*log(2.0)).*(-Ti-Tlead+t).^((Tc.*(-log(2.0)))./(Tc-Thr+Tc.*log(Thr./Tc)))+Fmax.*exp((Tc.*(log(Tc)-1.0).*log(2.0))./(Tc-Thr+Tc.*log(Thr./Tc))).*exp(((Ti+Tlead-t).*(-log(2.0)))./(Tc-Thr+Tc.*log(Thr./Tc))).*log(Thr./Tc).*(-Ti-Tlead+t).^((Tc.*(-log(2.0)))./(Tc-Thr+Tc.*log(Thr./Tc))).*(Ti+Tlead-t).*1.0./(Tc-Thr+Tc.*log(Thr./Tc)).^2.*log(2.0);
dFdThr=Fmax.*exp((Tc.*(log(Tc)-1.0).*log(2.0))./(Tc-Thr+Tc.*log(Thr./Tc))).*exp(((Ti+Tlead-t).*(-log(2.0)))./(Tc-Thr+Tc.*log(Thr./Tc))).*(Tc./Thr-1.0).*(-Ti-Tlead+t).^((Tc.*(-log(2.0)))./(Tc-Thr+Tc.*log(Thr./Tc))).*(Ti+Tlead-t).*1.0./(Tc-Thr+Tc.*log(Thr./Tc)).^2.*log(2.0)-Fmax.*Tc.*exp((Tc.*(log(Tc)-1.0).*log(2.0))./(Tc-Thr+Tc.*log(Thr./Tc))).*exp(((Ti+Tlead-t).*(-log(2.0)))./(Tc-Thr+Tc.*log(Thr./Tc))).*(log(Tc)-1.0).*(Tc./Thr-1.0).*(-Ti-Tlead+t).^((Tc.*(-log(2.0)))./(Tc-Thr+Tc.*log(Thr./Tc))).*1.0./(Tc-Thr+Tc.*log(Thr./Tc)).^2.*log(2.0)+Fmax.*Tc.*exp((Tc.*(log(Tc)-1.0).*log(2.0))./(Tc-Thr+Tc.*log(Thr./Tc))).*exp(((Ti+Tlead-t).*(-log(2.0)))./(Tc-Thr+Tc.*log(Thr./Tc))).*log(-Ti-Tlead+t).*(Tc./Thr-1.0).*(-Ti-Tlead+t).^((Tc.*(-log(2.0)))./(Tc-Thr+Tc.*log(Thr./Tc))).*1.0./(Tc-Thr+Tc.*log(Thr./Tc)).^2.*log(2.0);
dFdFmax=exp((Tc.*(log(Tc)-1.0).*log(2.0))./(Tc-Thr+Tc.*log(Thr./Tc))).*exp(((Ti+Tlead-t).*(-log(2.0)))./(Tc-Thr+Tc.*log(Thr./Tc))).*(-Ti-Tlead+t).^((Tc.*(-log(2.0)))./(Tc-Thr+Tc.*log(Thr./Tc)));

t=linspace(Ti+Tlead,Ti+Tlead+Ttot,Ttot*fs*10^(-3));