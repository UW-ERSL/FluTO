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
%% squircle
disp('generating squircle ...')
mstrIdentifier = 0;
ctr = 1; % keep track of how many have been added
numA = 100; % num samples about "a" length of rect
a = linspace(0.01, 0.99, numA); % 0.95 to avoid full blocking
for i = 1:numA
    
        [mstr, area_, length_] = squircle(nelx, nely, a(i),0.11);
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

% %% fish_scale 
% disp('generating fishScale ...')
% mstrIdentifier = mstrIdentifier + 1;
% numA = 100; % num samples about "a" length of rect
% a = linspace(0.01, 0.99, numA); % 0.95 to avoid full blocking
% for i = 1:numA
% 
%         [mstr, area_, length_] = fish_scale_super_ellipse(nelx, nely, a(i), 0.316);
% %         ch = fluidHomogenization(lx,ly,matProp,cellAngleDeg,mstr);
%         microstructures(mstrIdentifier+1,i,:,:) = mstr;
% %         identifier(ctr) = mstrIdentifier; % class of rect microstr
% %         mstrsize(ctr) = a(i); % append area to vec
% %         length(ctr) = length_; % append len to vec
% %         c00(ctr) = ch(1,1);
% %         c11(ctr) = ch(2,2);
% %         c01(ctr) = ch(1,2);
% %         ctr = ctr + 1;
%    
% end

%% fish_scale_half 
disp('generating fishScale half ratio ...')
mstrIdentifier = mstrIdentifier + 1;
numA = 100; % num samples about "a" length of rect
a = linspace(0.01, 0.99, numA); % 0.95 to avoid full blocking
for i = 1:numA

        [mstr, area_, length_] = fish_scale_super_ellipse_half(nelx, nely, a(i), 0.658);
%         ch = fluidHomogenization(lx,ly,matProp,cellAngleDeg,mstr);
        microstructures(mstrIdentifier+1,i,:,:) = mstr;
%         identifier(ctr) = mstrIdentifier; % class of rect microstr
%         mstrsize(ctr) = a(i); % append area to vec
%         length(ctr) = length_; % append len to vec
%         c00(ctr) = ch(1,1);
%         c11(ctr) = ch(2,2);
%         c01(ctr) = ch(1,2);
%         ctr = ctr + 1;
%    
end
% %% square
disp('generating square ...')
mstrIdentifier = mstrIdentifier + 1;
numA = 100; % num samples about "a" length of rect
a = linspace(0.01, 0.99, numA); % 0.95 to avoid full blocking
for i = 1:numA
 
        [mstr, area_, length_] = square(nelx, nely, a(i),0);
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

% % %% triangle
% disp('generating triangle ...')
% mstrIdentifier = mstrIdentifier + 1;
% numA = 100; % num samples about "a" length of rect
% a = linspace(0.01, 0.99, numA); % 0.95 to avoid full blocking
% for i = 1:numA
%  
%         [mstr, area_, length_] = Triangle(nelx, nely, a(i),0.567);
% %         ch = fluidHomogenization(lx,ly,matProp,cellAngleDeg,mstr);
%         microstructures(mstrIdentifier+1,i,:,:) = mstr;
% %         identifier(ctr) = mstrIdentifier; % class of rect microstr
% %         mstrsize(ctr) = a(i); % append area to vec
% %         length(ctr) = length_; % append len to vec
% %         c00(ctr) = ch(1,1);
% %         c11(ctr) = ch(2,2);
% %         c01(ctr) = ch(1,2);
% %         ctr = ctr + 1;
%   
% end


% %% circle
disp('generating circle ...')
mstrIdentifier = mstrIdentifier + 1;
numA = 100; % num samples about "a" length of rect
a = linspace(0.01, 0.99, numA); % 0.95 to avoid full blocking
for i = 1:numA
 
        [mstr, area_, length_] = circle(nelx, nely, a(i), 0.215);
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

% %% ellipse
disp('generating ellipse ...')
mstrIdentifier = mstrIdentifier + 1;
numA = 100; % num samples about "a" length of rect
a = linspace(0.01, 0.99, numA); % 0.95 to avoid full blocking
for i = 1:numA
 
        [mstr, area_, length_] = ellipse(nelx, nely, a(i), 0.607);
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

% % %% plus
% disp('generating plus ...')
% mstrIdentifier = mstrIdentifier + 1;
% numA = 100; % num samples about "a" length of rect
% a = linspace(0.01, 0.99, numA); % 0.95 to avoid full blocking
% for i = 1:numA
%  
%         [mstr, area_, length_] = Plus(nelx, nely, a(i), 4/9);
% %         ch = fluidHomogenization(lx,ly,matProp,cellAngleDeg,mstr);
%         microstructures(mstrIdentifier+1,i,:,:) = mstr;
% %         identifier(ctr) = mstrIdentifier; % class of rect microstr
% %         mstrsize(ctr) = a(i); % append area to vec
% %         length(ctr) = length_; % append len to vec
% %         c00(ctr) = ch(1,1);
% %         c11(ctr) = ch(2,2);
% %         c01(ctr) = ch(1,2);
% %         ctr = ctr + 1;
%   
% end

% %% mucosa 10
% disp('generating mucosa 10 ...')
% mstrIdentifier = mstrIdentifier + 1;
% numA = 100; % num samples about "a" length of rect
% a = linspace(0.01, 0.99, numA); % 0.95 to avoid full blocking
% for i = 1:numA
%  
%         [mstr, area_, length_] = mucosa(nelx, nely, a(i), 0.2901, 10);
% %         ch = fluidHomogenization(lx,ly,matProp,cellAngleDeg,mstr);
%         microstructures(mstrIdentifier+1,i,:,:) = mstr;
% %         identifier(ctr) = mstrIdentifier; % class of rect microstr
% %         mstrsize(ctr) = a(i); % append area to vec
% %         length(ctr) = length_; % append len to vec
% %         c00(ctr) = ch(1,1);
% %         c11(ctr) = ch(2,2);
% %         c01(ctr) = ch(1,2);
% %         ctr = ctr + 1;
%   
% end

% %% mucosa 20
disp('generating mucosa 20 ...')
mstrIdentifier = mstrIdentifier + 1;
numA = 100; % num samples about "a" length of rect
a = linspace(0.01, 0.99, numA); % 0.95 to avoid full blocking
for i = 1:numA
 
        [mstr, area_, length_] = mucosa(nelx, nely, a(i), 0.2919, 20);
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
