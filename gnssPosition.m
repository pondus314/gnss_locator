function [satPos] = gnssPosition(gnssSvTime,gnssSvWeek,gnssSvIds,gnssSvEph)
%gnssPosition Calculates the position of a gnss satellite in ECEF
%coordinates from its ephemeris at received time of transmission
%   inputs: gnssSvTime - received time in seconds
gpsEph = ClosestGpsEph(gnssSvEph,gnssSvIds,gnssSvTime(1)+gnssSvWeek(1)*GpsConstants.WEEKSEC);
dtsvS = gpsDtsv(gpsEph,gnssSvTime+gnssSvWeek*GpsConstants.WEEKSEC);
t_gps = gnssSvTime - dtsvS.';

[svXyzTtx,dtsv,svXyzDot,dtsvDot]=GpsEph2Pvt(gpsEph.',[gnssSvWeek,t_gps]);

satPos = [svXyzTtx, t_gps];
end
