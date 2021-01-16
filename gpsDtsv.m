function [dtsvS] = gpsDtsv(gpsEph,tS)
%GPSDTSV Summary of this function goes here
%   Detailed explanation goes here



%% Check size of inputs
if min(size(tS))>1
  error('tS must be a vector or a scalar, not a matrix')
end
tS=tS(:)';%make t a row vector, to match [gpsEph.*] vectors below
pt=length(tS);
[p]=length(gpsEph);
if (p>1 && pt~=p), 
  error('If gpsEph is a vector tS must be a vector with #rows = length(gpsEph),\n')
end
%%
%%Extract the necessary variables from gpsEph 
TGD     = [gpsEph.TGD];
Toc     = [gpsEph.Toc];
af2     = [gpsEph.af2];
af1     = [gpsEph.af1];
af0     = [gpsEph.af0];
Delta_n = [gpsEph.Delta_n];
M0      = [gpsEph.M0];
e       = [gpsEph.e];
Asqrt   = [gpsEph.Asqrt];
Toe     = [gpsEph.Toe];
%%
%%Calculate dependent variables ------------------------------------------------

tk = tS - Toe; %time since time of applicability
I = find(tk > 302400.0);
if any(I),  tk(I) = tk(I)-GpsConstants.WEEKSEC; end,
I = find(tk < -302400.0);
if (I), tk(I) = tk(I)+GpsConstants.WEEKSEC; end,

A = Asqrt.^2;    %semi-major axis of orbit
n0=sqrt(GpsConstants.mu./(A.^3));    %Computed mean motion (rad/sec)
n=n0+Delta_n;      %Corrected Mean Motion
Mk=M0+n.*tk;     %Mean Anomaly
Ek=Kepler(Mk,e);  %Solve Kepler's equation for eccentric anomaly
%%
%%Calculate satellite clock bias (See ICD-GPS-200 20.3.3.3.3.1) ----------------
dt = tS - Toc;
I = find(dt > 302400.0);
if any(I),  dt(I) = dt(I)-GpsConstants.WEEKSEC; end
I = find(dt < -302400.0);
if (I), dt(I) = dt(I)+GpsConstants.WEEKSEC; end

dtsvS = af0 + af1.*dt + af2.*(dt.^2)  + ...
    GpsConstants.FREL.*e.*Asqrt.*sin(Ek) -TGD;

end % end of function GpsEph2Dtsv


