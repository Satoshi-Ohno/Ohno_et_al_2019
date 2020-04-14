function  [randSeed] = init(inputRandSeed)
%% Random number generation

if nargin==0
    inputRandSeed = 1;
end
if isempty(inputRandSeed)
    inputRandSeed = 1;
end
if isnan(inputRandSeed)
    inputRandSeed = 1;    
end

currentTime = clock; 
second = currentTime(6); 

randSeed = second*inputRandSeed*1000; 
RandStream.setGlobalStream(RandStream('mt19937ar', 'Seed', randSeed)); 

%%  record variables except for randSeed
varList = [];
varList = who;
varList = setdiff(varList, {'randSeed'});

%% initialize
clc
clear(varList{:})
close all
tic
disp(datestr(clock))

end
