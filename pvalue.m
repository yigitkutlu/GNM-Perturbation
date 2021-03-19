clc; clear all;

%% INPUT RULES % FILES

% sections C alpha calculations & All atom calculations

% related region's pdb file (all atom)
protein = pdbread('HRasGTPase.pdb');

% related region's ca only file - to overcome missing residue issues
fname1='HRasGTPase_ca.pdb'; % alpha carbon coordinates of the molecule
[d1, d2, d3, d4, d5, d6, x_1, y_1, z_1, d7, beta, d8]=textread(fname1,'%s%f%s%s%s%f%f%f%f%f%f%s');

resnum=length(d6)-1;

k=1;
% minima obtained from GNM_peturbation
% "locs" will give you the minima
min = [4 6 10 16 19 24 40 42 46 53 55 60 68 71 76 78 81 89 93 96 100 109 111 114 137 139 142 152 157 160]; %minima positions

% related data (obtained from deep sequencing studies)
DATA = [6 10 15 16 17 35 38 40 57 75 78 111 112 113 114 115 119 138 141 156 159]; % Important positions from data

% IMPORTANT: this is important positions
% If real residue index is used change the following
% D=d6(DATA(i)); is not needed
% if (protein.Model.Atom(n).resSeq == D); change it to 
% if (protein.Model.Atom(n).resSeq == DATA(i));

%%  ----- CA ANALYSIS ----------

% for n=1:length(protein.Model.Atom);
%    if strcmp(protein.Model.Atom(n).AtomName,'CA');
%       CA(k,:) = [protein.Model.Atom(n).X protein.Model.Atom(n).Y protein.Model.Atom(n).Z];
%       k = k+1;
%    end    
% end
% 
% len=length(protein.Model.Atom);
%  
% for i=1:length(min)
%     A=min(i);
%     MinCA(i,:) = CA(A,:);
% end
%     
% Randist = ones(size(DATA))*100;
% 
% for n=1:10000
%     r1 = randperm(resnum,length(DATA));
%     Randist1(n,:) = Randist;
%     for i=1:length(DATA)
%     B=r1(i);
%     RanCA(i,:) = CA(B,:);
%       for j=1:length(min)
%        Dist = sqrt(((RanCA(i,1))-(MinCA(j,1)))^2 + ((RanCA(i,2))-(MinCA(j,2)))^2 + ((RanCA(i,3))-(MinCA(j,3)))^2);
%         if Dist < Randist1 (n,i)
%            Randist1 (n,i) = Dist;
%         end
%       end
%     end
%   totaldist(n) = sum(Randist1 (n,:));    
% end
% 
% meandist = totaldist/length(DATA);
% 
% 
% Resdist = ones(size(DATA))*100;
%  for i=1:length(DATA)
%     res=DATA(i);
%     resCA(i,:) = CA(res,:);
%       for j=1:length(min)
%        Dist = sqrt(((resCA(i,1))-(MinCA(j,1)))^2 + ((resCA(i,2))-(MinCA(j,2)))^2 + ((resCA(i,3))-(MinCA(j,3)))^2);
%         if Dist < Resdist (i)
%            Resdist (i) = Dist;
%         end
%       end
%  end
% 
% restotaldist = sum(Resdist);
% resmeandist = restotaldist/length(DATA);
%  
% Zer1 = Randist1';
% zeros = sum(Zer1 ==0);
% 
% Zerres = Resdist';
% reszeros = sum(Zerres ==0);
% DATAz=reszeros; % number of overlapping positions to minima with important positions

%% plot the histogram

% figure (1)
% h1=histfit(zeros,13);

%% calculation for overlapping positions  

% m1= mean(zeros); %mean number of overlapping positions to minima with random positions
% s1= std(zeros);
% z1= ((DATAz-m1)/s1);
% p1 = 1-tcdf(z1,9999);

%% calculation for C-alpha distance 

% m2 = mean(totaldist); %Total mean distance of minima to random positions (C alpha)
% s2 = std(totaldist);
% z2= ((m2-restotaldist)/s2);
% p2 = 1-tcdf(z2,9999);
% 
% m3 = mean(meandist); % mean average distance of minima to random positions (C alpha)
% s3 = std(meandist);
% z3= ((m3-resmeandist)/s3);
% p3 = 1-tcdf(z3,9999);
% figure (2)
% h2= histfit(meandist,13);
% figure (3)
% h3= histfit(totaldist,13);


%% ---- ALL ATOM ANALYSIS ---- 

for i=1:length(min)
    for n=1:length(protein.Model.Atom);
      if (protein.Model.Atom(n).resSeq == d6((min(i))));
         minres(k,:) = [protein.Model.Atom(n).X protein.Model.Atom(n).Y protein.Model.Atom(n).Z];
         k = k+1;
      end    
    end
end

Randistres = ones(size(DATA))*100;
for f=1:10000
     n=1;
     r2 = randperm(resnum,length(DATA));
     Randistres1(f,:) = Randistres;
     for i=1:length(DATA)
         k=1;
         randres=zeros(k,3);
         C=d6(r2(i));
         for n=1:length(protein.Model.Atom);
            if (protein.Model.Atom(n).resSeq == C);
            randres(k,:) = [protein.Model.Atom(n).X protein.Model.Atom(n).Y protein.Model.Atom(n).Z];
            k = k+1;
            end    
         end
         for m=1:length(randres)
             for j=1:length(minres)
              Dist = sqrt(((randres(m,1))-(minres(j,1)))^2 + ((randres(m,2))-(minres(j,2)))^2 + ((randres(m,3))-(minres(j,3)))^2);
               if Dist < Randistres1 (f,i)
                Randistres1 (f,i) = Dist;
               end
             end
         end
     end
  totaldistres(f) = sum(Randistres1 (f,:));    
end

Selresdist = ones(size(DATA))*100; %large numbers for distance

for i=1:length(DATA)
    k=1;
    selres=zeros(k,3);
    D=d6(DATA(i));
    for n=1:length(protein.Model.Atom);
      if (protein.Model.Atom(n).resSeq == D);
         selres(k,:) = [protein.Model.Atom(n).X protein.Model.Atom(n).Y protein.Model.Atom(n).Z];
         k = k+1;
      end    
    end
   for m=1:length(selres)
       for j=1:length(minres)
           Dist = sqrt(((selres(m,1))-(minres(j,1)))^2 + ((selres(m,2))-(minres(j,2)))^2 + ((selres(m,3))-(minres(j,3)))^2);
            if Dist < Selresdist (i)
                Selresdist (i) = Dist;
            end
       end
   end
end

meandistres=totaldistres/length(DATA);

seltotaldist= sum(Selresdist); %total distance of minima to important positions
selmeandist= seltotaldist/length(DATA);  %mean distance of minima to important positions

%% p-value calculations

%p-value for total distance
m4 = mean(totaldistres); % total mean distance of minima to random positions
s4 = std(totaldistres);
z4= ((m4-seltotaldist)/s4);
p4 = 1-tcdf(z4,9999);

%p-value for mean distance (p-values should be equal)
m5 = mean(meandistres);  %mean average distance of minima to random positions
s5 = std(meandistres);
z5= ((m5-selmeandist)/s5);
p5 = 1-tcdf(z5,9999);

%% plot the histogram
figure(2)
h4= histfit(meandistres,15);