clear all; close all;clc;
%% Parameters
nelx = 100;
nely = 100;
lx = 1; % length of microstr along x
ly = 1;
solidPermeability = 1e6;
fluidPermeability = 0;
matProp = [solidPermeability fluidPermeability];
cellAngleDeg = 90;
%% initalization
microstructures = [];
mstrsize = [];
length = [];
c00 = [];
c11 = [];
c10 = [];
identifier = [];

% %% mucosa 20
disp('generating mucosa 20 ...')
mstrIdentifier = 0;
numA = 100; % num samples about "a" length of rect
a = linspace(0.01, 0.99, numA); % 0.95 to avoid full blocking
for i = 1:numA
 
        [mstr, area_, length_] = mucosa(nelx, nely, a(i), 0.2901, 10);
%         ch = fluidHomogenization(lx,ly,matProp,cellAngleDeg,mstr);
        microstructures(mstrIdentifier+1,i,:,:) = mstr;
%         identifier(ctr) = mstrIdentifier; % class of rect microstr
%         mstrsize(ctr) = a(i); % append area to vec
%         length(ctr) = length_; % append len to vec
%         c00(ctr) = ch(1,1);
%         c11(ctr) = ch(2,2);
%         c01(ctr) = ch(1,2);
%         ctr = ctr + 1;
  
end

% %% void
disp('generating void ...')
mstrIdentifier = mstrIdentifier + 1;
numA = 100; % num samples about "a" length of rect
a = linspace(0.01, 0.99, numA); % 0.95 to avoid full blocking
for i = 1:numA
 
        [mstr, area_, length_] = square(nelx, nely, a(i),0.9999);
%         ch = fluidHomogenization(lx,ly,matProp,cellAngleDeg,mstr);
        microstructures(mstrIdentifier+1,i,:,:) = mstr;
%         identifier(ctr) = mstrIdentifier; % class of rect microstr
%         mstrsize(ctr) = a(i); % append area to vec
%         length(ctr) = length_; % append len to vec
%         c00(ctr) = ch(1,1);
%         c11(ctr) = ch(2,2);
%         c01(ctr) = ch(1,2);
%         ctr = ctr + 1;
  
end

%% write to file

% T = table(identifier', mstrsize',length', c00', c11');
% filename = 'fluidMicrostructureData.xlsx';
% writetable(T,filename,'Sheet',1,'WriteVariableNames',true);
% 
