function [pulse,azerr,elerr] = getNearestUCDpulse(azimuth,elevation,h3D)
%   retrieves the impulse response from h3D that is closest to the specified 
%   azimuth and elevation (in degrees)
%   h3D is the array containing all HRTFs for a given subject (left or right)


elmax = 50;
elindices = 1:elmax;
elevations = -45 + 5.625*(elindices-1);
azimuths = [-80 -65 -55 -45:5:45 55 65 80];

%el est est l'indice du plus proche de l'indicie elevation entré
el = round ((elevation+45)/5.625 +1);
el= max(el,1);
el = min (el,elmax);
elerr= el - (elevation+45)/5.625 +1;


%azim est l'indice du plus proche de l'indice azimuth entré
[azerr, azim] = min (abs(azimuths - azimuth));


fprintf('les paramètres utilisés pour la HRTF sont:\n \t azimuth: %i degrés \n \t élevation: %6.3f degrés',azimuths(azim),-45 + 5.625*(el-1));
fprintf('\nles indices sont:\n \t azim: %i \n \t el: %d ',azim,el);
azimuths(azim);

pulse = squeeze(h3D(azim,el,:));
