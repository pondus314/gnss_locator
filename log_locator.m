%% Raw data import
opts = delimitedTextImportOptions;
opts.DataLines = 12;
opts.VariableNamesLine = 6;
filein = readcell("gnss_log_2020_12_17_12_18_20.txt",opts);

raw_rows = any(strcmp(filein, 'Raw'), 2);
raws_unfilt = filein(raw_rows, 2:31);

for i  = [1:3]
    raws_unfilt(:,23) = [];
end

raws_unfilt(:,[3,4,24,27]) = {0}; % used when reading files from the new google logger

raws_mat = cell2mat(raws_unfilt);
%% Raw data filtering and grouping
state_filter = (bitand(raws_mat(:,13),1) & bitand(raws_mat(:,13),8));
constell_filter = (raws_mat(:,25) == 1);
sv_uncert_filter = (raws_mat(:,15) < 500);
raw_filter = state_filter & constell_filter & sv_uncert_filter;
raws = raws_mat(raw_filter,:);

[satellite,svids,constells] = findgroups(raws(:,11),raws(:,25));
[measurement, times] = findgroups(-raws(:,5));
N = numel(times)
M = numel(svids)
K = size(raws,1)
L = size(raws,2)
%% Satellite position calculation
utcTime = [2020, 12, 17, 12, 18, 25];
allGpsEph = GetNasaHourlyEphemeris(utcTime,'eph');

gpsFilter = constells(satellite)==1;
% remove the +1 in the next line if using the google logger
satPos = gnssPosition(raws(gpsFilter,14)*10.^-9,-floor(raws(gpsFilter,5)/GpsConstants.WEEKSEC*10.^-9), raws(gpsFilter, 11),allGpsEph)
%% Measurement organisation
meas_organised = zeros(N,M,size(raws,2)+4);
measured = false(N,M);

for k = [1:K]
    n = measurement(k);
    m = satellite(k);
    meas_organised(n, m, 1:L) = raws(k, :);
    for i  = [1:4]
        meas_organised(n,m,L+i) = satPos(k,i);
    end
    measured(n, m) = true;
end

jtNs = 2;            % time Nanos
jtFBNs = 5;          % Full Bias Nanos 
jtBns = 6;           % Bias Nanos
jtONs = 12;          % Time offset Nanos
jsvtUNs = 15;        % rtx uncertainty Nanos
jsvXYZM = (L+1:L+3); % sv ECEF Position XYZ [m]
jsvTS  = L+4;        % sv Gps Time [s]
"done"
%% User position calculation
userPoss = zeros(4,N);
for i = [1:N]
    meas_i = reshape(meas_organised(i, measured(i,:),:), [], L+4);
    svPos=meas_i(:,[jsvXYZM,jsvTS]);
    if numel(svPos) < 16
        continue
    end
    userTime = mod((meas_i(1, jtFBNs)+meas_i(1, jtNs))*10.^-9,GpsConstants.WEEKSEC);
    userPos = userPosition(svPos(1:4,:),userTime);
    userPoss(:,i) = userPos;
end

pos =ecef2lla(userPoss(1:3,:).');

pos(:,1:2)

%% backwards check
userecef = [3916732, 8346, 5017023 ];
dif = svPos(:,1:3)-userecef;

geometric_r = vecnorm(dif,2,2);

pr = (userTime - svPos(:,4)) * GpsConstants.LIGHTSPEED; % 
pr - geometric_r % difference between the measured range and expected range in metres

%% checking that the equation solution works
meas = reshape(meas_organised(1, measured(1,:),:), [], L+4);
svPos=meas(:,[jsvXYZM,jsvTS]);
userTime = mod((-meas(1,jtFBNs)+meas(1,jtNs))*10.^-9,GpsConstants.WEEKSEC);
userPos = userPosition(svPos(1:4,:),userTime)
expPos = lla2ecef([52.211094, 0.091276, 10])

realTime = userPos(4) % user time found by solving the gps equations
distmat = svPos(1:4,[1:3])-userPos([1:3]) % distance matrix from individual satellites
grs = vecnorm(distmat,2,2)/GpsConstants.LIGHTSPEED % geometric range to satellites in seconds
rxExpTS = svPos(1:4,4)+grs % expected time of signal reception
(rxExpTS -  realTime) % should be four zeros
