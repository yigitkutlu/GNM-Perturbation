clc; clear all;

%% type in the input file name/location

% read single pdb file which contains alpha-carbons only!!!!
% the file should not include any TER or ANISOU!!!
% only only one atomic coordinate per amino acid (C alpha, C beta etc.)

% for mode selection see "mode selection" section
% for fiugre adjustments see plot sections

fname1='pdz_ca.pdb'; % alpha carbon coordinates of the molecule
[d1, d2, d3, d4, d5, d6, x_1, y_1, z_1, d7, beta, d8]=textread(fname1,'%s%f%s%s%s%f%f%f%f%f%f%s');

%% 
x=[x_1']; 
y=[y_1'];
z=[z_1'];
beta_fact=[beta];

% erase the last line
x=[x(1:size(x,2)-1)];
y=[y(1:size(y,2)-1)];
z=[z(1:size(z,2)-1)];
beta=[beta(1:size(beta,1)-1)];

clear *_*;
clear d1 d2 d3 d4 d5 d6 d7 d8;

%  general
resnum=size(x,2);
frame=size(x,1);

rcut_gnm=10;

% gamma constant
ga=1;

%% mode selection
% select the number of modes ie. 3, 5, 10 etc. for all modes use "resnum-1"

modes = 3; 
%modes = resnum -1; % all modes

%% preallocate arrays
r2=zeros(resnum,resnum);
cont=zeros(resnum,resnum);
invcont=zeros(resnum,resnum);
crosscorr=zeros(resnum,resnum);
invcont1=zeros(resnum,resnum);
crosscorr1=zeros(resnum,resnum);
invcont_av10=zeros(resnum,resnum);
crosscorr_av10=zeros(resnum,resnum);

%diagonal=zeros(resnum);
U=zeros(resnum,resnum);
S=zeros(resnum,resnum);
V=zeros(resnum,resnum);
w_snsh=zeros(resnum,frame);
w1_snsh=zeros(resnum,frame);
w2_snsh=zeros(resnum,frame);
w3_snsh=zeros(resnum,frame);



%% Unperturbed calculations
for snapshot=1:frame
      
        for j=1:resnum
            for k=1:resnum
                distx = x(snapshot,j)-x(snapshot,k);   
                disty = y(snapshot,j)-y(snapshot,k);   
                distz = z(snapshot,j)-z(snapshot,k);
                
                r2(j,k)=distx^2+disty^2+distz^2;
                
                r=sqrt(distx^2+disty^2+distz^2);
                
                % Kirchhoff-Connectivity matrix
              
                if (r <= rcut_gnm && j~=k && r > 0.0001)
                    
                    cont(j,k)=-1;
                    
                else
                    cont(j,k)=0;                                
                end             
            end
        end
       
        % detailed balance for connectivity matrix
        diagonal=sum(cont(:,:));
        
        for j=1:resnum
            for i=1:resnum
                if i == j
                    cont(i,j)=-1*diagonal(i); % Kirchhoff
                end
            end
        end
    
                
        [U(:,:),S(:,:),V(:,:)]=svd(cont(:,:));
        w=diag(S(:,:));
        eign1=w;
        
        % to export as output
        w_snsh(:,snapshot)=diag(S(:,:));        
        
	    w1_snsh(snapshot)=w_snsh(resnum-1,snapshot);
        w2_snsh(snapshot)=w_snsh(resnum-2,snapshot);
        w3_snsh(snapshot)=w_snsh(resnum-3,snapshot);
		   
%         % GNM all eigenvectors and modes

        for j=resnum-1:-1:resnum-modes %slowest n modes
            for i=1:resnum
                slow_vectors_gnm(i,j)=V(i,j);
                slow_modes_gnm(i,j)=V(i,j)*V(i,j)/w(j);
            end
        end
        
        
        % GNM all eigenvectors and modes overall contribution
        for i=1:resnum %for all modes
            
            cummulative_slow_vec(i)=sum(slow_vectors_gnm(i,:));
            cummulative_slow_modes(i)=sum(slow_modes_gnm(i,:));
            
        end
        cummulative_slow_vec=cummulative_slow_vec';
        cummulative_slow_modes=cummulative_slow_modes';
        
end

%% Perturbed Calculations

pert=cont;
p=1;
t=2;

for c=1:resnum

    for k=1:resnum
        
        if (pert(c,k)==-1)
            pert(c,k)=-1.2;
        end        
    end
    for j=1:resnum
        for k=1:resnum
            if(j==k)
                pert(j,k)=0;
            end
        end
    end
    
    % detailed balance for connectivity matrix
    diagonal=sum(pert(:,:));
        
    for j=1:resnum
        for i=1:resnum
        	if (i == j)
                pert(i,j)=-1*diagonal(i); % Kirchhoff
            end
        end
    end
    
    [U(:,:),S(:,:),V(:,:)]=svd(pert(:,:));
    w=diag(S(:,:));
    eign(:,c)=w;
		   
   
     for j=resnum-1:-1:resnum-modes  % slowest n modes
        for i=1:resnum
            prt_slow_vectors_gnm(i,j)=V(i,j);
            prt_slow_modes_gnm(i,j)=V(i,j)*V(i,j)/w(j);
        end
    end
    
    
    for i=1:resnum   % GNM all eigenvectors and modes overall contribution
        prt_cummulative_slow_vec(i)=sum(prt_slow_vectors_gnm(i,:));
        prt_cummulative_slow_modes(i)=sum(prt_slow_modes_gnm(i,:));          
    end
    prt_cummulative_slow_vec1=prt_cummulative_slow_vec.';
    prt_cummulative_slow_modes1=prt_cummulative_slow_modes.';
    
    p=p+1;

    
    %Calculates the differences in fluctuations after perturbation and
    %creates the modes_diff matrix (first column: residue index)
    for i=1:resnum
     
           D(i)=(prt_cummulative_slow_modes1(i)-cummulative_slow_modes(i)); %fluctuation difference
           
    end
    
    modes_diff(1:resnum,c)=D;
    
    t=t+1;
    scan1(:,:,c)=pert(:,:);
    pert=cont; 
end

tot_diff= sum(modes_diff); %cumulative fluctuation difference of each perturbed residue

%% eigenvalue calculations

for i=resnum-1:-1:resnum-modes
    
    eign2(resnum-i,1)=eign1(i);
    for j=1:resnum
        eignpert(resnum-i,j) = eign(i,j);
    end
end

pertdif = eignpert-eign2;
tpertdiff= sum(pertdif);  %cumulative eigenvalue difference of each perturbed residue

%% cumulative fluctuation difference plot

% in order to mark important residues on the figures

%important positions

X1=[13 14 15 17 18 19 20 26 28 31 37 43 49 52 57 62 65 66 69 78];

%Real residue index of important positions

str1={' 323', ' 324', ' 325', ' 327', ' 328', ' 329', ' 330', ' 336', ' 338', ' 341', ' 347', ' 353', ' 359', ' 362', ' 367', ' 372', ' 375', ' 376', ' 379', ' 388'}; %important res

% minor important positions (binding site etc.

 X2=[ 12 13 14 15 16 17 18 29 62 69 70]; %binding positions
% X3=[ 12 13 14 15 16 17 18 29 62 69 70]; %other positions

figure(1)
axes1 = axes;
box(axes1,'on');
set(axes1,'FontSize',22);
x0=10;
y0=10;
width=1920;
height=963;
set(gcf,'position',[x0,y0,width,height])
hold(axes1,'on');
plot(tot_diff,'DisplayName','Residues',...
    'MarkerFaceColor',[0 0.447058826684952 0.74117648601532],...
    'MarkerEdgeColor',[0 0.447058826684952 0.74117648601532],...
    'MarkerSize',4,...
    'Marker','hexagram',...
    'LineWidth',1);

% type in area of interest.

xlim([0 89])
%ylim([-0.03 0.025])

% change title or labels as pleased.

xlabel('Perturbed Residue Index','FontSize',26,'FontWeight','bold')
ylabel('Cumulative Fluctuation Difference','FontSize',26,'FontWeight','bold')
title('Three slowest modes of the PSD95-PDZ domain','FontSize',30,'FontWeight','bold')

% set positions as real residue index

set(axes1,'XTick',[1 10 20 30 40 50 60 70 80],'XTickLabel',...
    {'311','320','330','340','350','360','370','380','390'});
hold on;

% plot important positions on the graph
%type in display names for correct labeling on the legend.

plot(X1,tot_diff(X1),'p','MarkerFaceColor','red','MarkerEdgeColor','red','MarkerSize',16,'DisplayName','Mutation Sensitive Positions');
text(X1,tot_diff(X1),str1,'FontSize',22,'Color','red','FontWeight','bold');
hold on 

plot(X2,tot_diff(X2),'o','MarkerFaceColor','g','MarkerEdgeColor','green','MarkerSize',6,'DisplayName','Known Binding Sites');

%plot(X3,tot_diff(X2),'o','MarkerFaceColor','g','MarkerEdgeColor','green','MarkerSize',6,'DisplayName','...');

% legend adjustments

legend1 = legend(axes1,'show');
set(legend1,...
    'Location','northeast',...
    'FontWeight','bold',...
    'FontSize',24);
hold off;

%% cumulative eigenvalue difference plot

% figure(2)
% axes1 = axes;
% box(axes1,'on');
% set(axes1,'FontSize',28);
% x0=10;
% y0=10;
% width=1920;
% height=963;
% set(gcf,'position',[x0,y0,width,height])
% hold(axes1,'on');
% plot(tpertdiff,'DisplayName','Residues',...
%     'MarkerFaceColor',[0 0.447058826684952 0.74117648601532],...
%     'MarkerEdgeColor',[0 0.447058826684952 0.74117648601532],...
%     'MarkerSize',4,...
%     'Marker','hexagram',...
%     'LineWidth',1);
% 
% % type in area of interest
% 
% xlim([0 89])
% %ylim([-0.03 0.025])
% 
% xlabel('Perturbed Residue Index','FontSize',32,'FontWeight','bold')
% ylabel('Cumulative Fluctuation Difference','FontSize',32,'FontWeight','bold')
% title('Three slowest modes of the PSD95-PDZ domain','FontSize',36,'FontWeight','bold')
% 
% %set positions as real residue index
% set(axes1,'XTick',[1 10 20 30 40 50 60 70 80],'XTickLabel',...
%     {'311','320','330','340','350','360','370','380','390'});
% hold on;
% 
% %plot important positions on the graph
% %type in display names for correct labeling on the legend.
% 
% plot(X1,tpertdiff(X1),'p','MarkerFaceColor','red','MarkerEdgeColor','red','MarkerSize',16,'DisplayName','Mutation Sensitive Positions');
% text(X1,tpertdiff(X1),str1,'FontSize',28,'Color','red','FontWeight','bold');
% hold on 
% 
% plot(X2,tpertdiff(X2),'o','MarkerFaceColor','g','MarkerEdgeColor','green','MarkerSize',6,'DisplayName','Known Binding Sites');
% 
% %plot(X3,tot_diff(X2),'o','MarkerFaceColor','g','MarkerEdgeColor','green','MarkerSize',6,'DisplayName','...');
% 
% 
% 
% legend1 = legend(axes1,'show');
% set(legend1,...
%     'Position',[0.718923614479392 0.828136865342873 0.158333329406257 0.0545171324709123],...
%     'FontWeight','bold',...
%     'FontSize',28);
% hold off;

%% minimum position calculations

tot_diff2=tot_diff*(-1000);
mean_diff = mean(tot_diff2);
[pks,locs] = findpeaks(tot_diff2,'MinPeakHeight',mean_diff);

% if needs adjustments -like specific domain- make the adjustments below

%example for PAB1 RRM2 domain specific region RRM2 domain

% for i=1:74
%        tot_diffa(i)= tot_diff(89+i);
% end
% 
% tot_diffa=tot_diffa*(-1000);
% mean_diff = mean(tot_diffa);
% [pks,locs] = findpeaks(tot_diffa,'MinPeakHeight',mean_diff);



 
