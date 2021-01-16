function [userPos, guesses] = userPosition(satPoss, userTimeGuess, userPosGuess)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if ~exist('userPosGuess', 'var') 
   userPosGuess = [0,0,0,0]; 
else 
    userPosGuess= [userPosGuess,0];
end
guesses = zeros(20,4);
N = size(satPoss,1);
prM = (userTimeGuess-satPoss(:,4))*GpsConstants.LIGHTSPEED;
for i = [1:20]
    rHatVec =  (satPoss(:,[1:3])-userPosGuess(1:3));
    rHat = vecnorm(rHatVec,2,2);
    H = [rHatVec./rHat,ones(N,1)];
    delRho = rHat - prM;
    guessDelta = H \ delRho;
    userPosGuess(1:3) = userPosGuess(1:3) + guessDelta(1:3).';
    userPosGuess(4) = guessDelta(4)/GpsConstants.LIGHTSPEED;
    guesses(i,:) = userPosGuess;
    if vecnorm(guessDelta(1:3)) < 0.01
        break;
    end
    
end
userPosGuess(4) = userPosGuess(4)+userTimeGuess;
guesses(:,4) = guesses(:,4);
userPos = userPosGuess;

end
