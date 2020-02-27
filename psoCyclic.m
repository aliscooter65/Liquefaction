function [ PSOResults ] = psoCyclic(Net)
disp('Particle Swarm Optimization Algorithm');
disp('Programmed By: Ali Hendi 2016');
disp(' ');
S0cy=14.13:(400-14.13)/9:400;



%%====== input data====
npop=160; %number of particles
nvar=5; %number of variables
maxit=20; %number of iterations
%% S0
TT=zeros(numel(S0cy),nvar);

for oo=1:numel(S0cy)
%%========================================
w=1;
wdamp=0.99;
c1=2;
c2=2;
%========s0 mean  (kPa)=======
xmin5=S0cy(oo);%14.13;
xmax5=S0cy(oo);%400;
dx5=xmax5-xmin5;
vmax5=0.1*dx5;
%======== Dr (%) =======
xmin6=0; 
xmax6=100;
dx6=xmax6-xmin6;
vmax6=0.1*dx6;
%======== FC (%) =======
xmin7=0; 
xmax7=100;
dx7=xmax7-xmin7;
vmax7=0.1*dx7;
%======== Cu =======
xmin8=1.52;
xmax8=28.12;
dx8=xmax8-xmin8;
vmax8=0.1*dx8;
%======== D50 (mm) =======
xmin9=0.029; 
xmax9=0.46;
dx9=xmax9-xmin9;
vmax9=0.1*dx9;
%============================cumulate all of restrictions ==============
xmin=[xmin5;xmin6;xmin7;xmin8;xmin9];
xmax=[xmax5;xmax6;xmax7;xmax8;xmax9];
vmax=[vmax5;vmax6;vmax7;vmax8;vmax9];
xmin=xmin';
xmax=xmax';
vmax=vmax';

empty_particle.position=[];
empty_particle.velocity=[];
empty_particle.cost=[];
empty_particle.pbest=[];
empty_particle.pbestcost=[];

particle=repmat(empty_particle,npop,1);

gbest=zeros(maxit,nvar);
gbestcost=zeros(maxit,1);

for it=1:maxit
    if it==1
        gbestcost(1)=inf;
        for i=1:npop
            particle(i).velocity=zeros(1,nvar);
            particle(i).position=xmin+(xmax-xmin).*rand(1,nvar);
            particle(i).cost=Cost1(particle(i).position,Net);
            particle(i).pbest=particle(i).position;
            particle(i).pbestcost=particle(i).cost;
            
            if particle(i).pbestcost<gbestcost(it)
                gbest(it,:)=particle(i).pbest;
                gbestcost(it)=particle(i).pbestcost;
            end
        end
    else
        gbest(it,:)=gbest(it-1,:);
        gbestcost(it)=gbestcost(it-1);
        for i=1:npop
            particle(i).velocity=w*particle(i).velocity...
                                +c1*rand*(particle(i).pbest-particle(i).position)...
                                +c2*rand*(gbest(it,:)-particle(i).position);
                            
            particle(i).velocity=min(max(particle(i).velocity,-vmax),vmax);
            
            particle(i).position=particle(i).position+particle(i).velocity;
            
            particle(i).position=min(max(particle(i).position,xmin),xmax);
            
            particle(i).cost=Cost1(particle(i).position,Net);
            
            if particle(i).cost<particle(i).pbestcost
                particle(i).pbest=particle(i).position;
                particle(i).pbestcost=particle(i).cost;

                if particle(i).pbestcost<gbestcost(it)
                    gbest(it,:)=particle(i).pbest;
                    gbestcost(it)=particle(i).pbestcost;
                end
            end
        end
    end
    
    
    w=w*wdamp;
    
end
disp(['Iteration min S0' num2str(oo) '/' num2str(numel(S0cy))]);
TT(oo,:)=gbest(end,:);
PSOResultsmin.S0.result=TT;
PSOResultsmin.S0.iter{oo}=gbestcost;
end
% %% Dr
% TT=zeros(numel(Drcy),nvar);
% 
% for oo=1:numel(Drcy)
% %%========================================
% w=1;
% wdamp=0.99;
% c1=2;
% c2=2;
% %========s0 mean  (kPa)=======
% xmin5=14.13;
% xmax5=400;
% dx5=xmax5-xmin5;
% vmax5=0.1*dx5;
% %======== Dr (%) =======
% xmin6=Drcy(oo);%0; 
% xmax6=Drcy(oo);%100;
% dx6=xmax6-xmin6;
% vmax6=0.1*dx6;
% %======== FC (%) =======
% xmin7=0; 
% xmax7=100;
% dx7=xmax7-xmin7;
% vmax7=0.1*dx7;
% %======== Cu =======
% xmin8=1.52;
% xmax8=28.12;
% dx8=xmax8-xmin8;
% vmax8=0.1*dx8;
% %======== D50 (mm) =======
% xmin9=0.029; 
% xmax9=0.46;
% dx9=xmax9-xmin9;
% vmax9=0.1*dx9;
% %============================cumulate all of restrictions ==============
% xmin=[xmin5;xmin6;xmin7;xmin8;xmin9];
% xmax=[xmax5;xmax6;xmax7;xmax8;xmax9];
% vmax=[vmax5;vmax6;vmax7;vmax8;vmax9];
% xmin=xmin';
% xmax=xmax';
% vmax=vmax';
% 
% empty_particle.position=[];
% empty_particle.velocity=[];
% empty_particle.cost=[];
% empty_particle.pbest=[];
% empty_particle.pbestcost=[];
% 
% particle=repmat(empty_particle,npop,1);
% 
% gbest=zeros(maxit,nvar);
% gbestcost=zeros(maxit,1);
% 
% for it=1:maxit
%     if it==1
%         gbestcost(1)=inf;
%         for i=1:npop
%             particle(i).velocity=zeros(1,nvar);
%             particle(i).position=xmin+(xmax-xmin).*rand(1,nvar);
%             particle(i).cost=Cost1(particle(i).position,Net);
%             particle(i).pbest=particle(i).position;
%             particle(i).pbestcost=particle(i).cost;
%             
%             if particle(i).pbestcost<gbestcost(it)
%                 gbest(it,:)=particle(i).pbest;
%                 gbestcost(it)=particle(i).pbestcost;
%             end
%         end
%     else
%         gbest(it,:)=gbest(it-1,:);
%         gbestcost(it)=gbestcost(it-1);
%         for i=1:npop
%             particle(i).velocity=w*particle(i).velocity...
%                                 +c1*rand*(particle(i).pbest-particle(i).position)...
%                                 +c2*rand*(gbest(it,:)-particle(i).position);
%                             
%             particle(i).velocity=min(max(particle(i).velocity,-vmax),vmax);
%             
%             particle(i).position=particle(i).position+particle(i).velocity;
%             
%             particle(i).position=min(max(particle(i).position,xmin),xmax);
%             
%             particle(i).cost=Cost1(particle(i).position,Net);
%             
%             if particle(i).cost<particle(i).pbestcost
%                 particle(i).pbest=particle(i).position;
%                 particle(i).pbestcost=particle(i).cost;
% 
%                 if particle(i).pbestcost<gbestcost(it)
%                     gbest(it,:)=particle(i).pbest;
%                     gbestcost(it)=particle(i).pbestcost;
%                 end
%             end
%         end
%     end
%     
%     
%     
%     w=w*wdamp;
%     
% end
% disp(['Iteration min Dr' num2str(oo) '/' num2str(numel(Drcy))]);
% TT(oo,:)=gbest(end,:);
% PSOResultsmin.Dr.result=TT;
% PSOResultsmin.Dr.iter{oo}=gbestcost;
% end
% %% FC
% TT=zeros(numel(Fccy),nvar);
% 
% for oo=1:numel(Fccy)
% %%========================================
% w=1;
% wdamp=0.99;
% c1=2;
% c2=2;
% %========s0 mean  (kPa)=======
% xmin5=14.13;
% xmax5=400;
% dx5=xmax5-xmin5;
% vmax5=0.1*dx5;
% %======== Dr (%) =======
% xmin6=0; 
% xmax6=100;
% dx6=xmax6-xmin6;
% vmax6=0.1*dx6;
% %======== FC (%) =======
% xmin7=Fccy(oo);%0; 
% xmax7=Fccy(oo);%100;
% dx7=xmax7-xmin7;
% vmax7=0.1*dx7;
% %======== Cu =======
% xmin8=1.52;
% xmax8=28.12;
% dx8=xmax8-xmin8;
% vmax8=0.1*dx8;
% %======== D50 (mm) =======
% xmin9=0.029; 
% xmax9=0.46;
% dx9=xmax9-xmin9;
% vmax9=0.1*dx9;
% %============================cumulate all of restrictions ==============
% xmin=[xmin5;xmin6;xmin7;xmin8;xmin9];
% xmax=[xmax5;xmax6;xmax7;xmax8;xmax9];
% vmax=[vmax5;vmax6;vmax7;vmax8;vmax9];
% xmin=xmin';
% xmax=xmax';
% vmax=vmax';
% 
% empty_particle.position=[];
% empty_particle.velocity=[];
% empty_particle.cost=[];
% empty_particle.pbest=[];
% empty_particle.pbestcost=[];
% 
% particle=repmat(empty_particle,npop,1);
% 
% gbest=zeros(maxit,nvar);
% gbestcost=zeros(maxit,1);
% 
% for it=1:maxit
%     if it==1
%         gbestcost(1)=inf;
%         for i=1:npop
%             particle(i).velocity=zeros(1,nvar);
%             particle(i).position=xmin+(xmax-xmin).*rand(1,nvar);
%             particle(i).cost=Cost1(particle(i).position,Net);
%             particle(i).pbest=particle(i).position;
%             particle(i).pbestcost=particle(i).cost;
%             
%             if particle(i).pbestcost<gbestcost(it)
%                 gbest(it,:)=particle(i).pbest;
%                 gbestcost(it)=particle(i).pbestcost;
%             end
%         end
%     else
%         gbest(it,:)=gbest(it-1,:);
%         gbestcost(it)=gbestcost(it-1);
%         for i=1:npop
%             particle(i).velocity=w*particle(i).velocity...
%                                 +c1*rand*(particle(i).pbest-particle(i).position)...
%                                 +c2*rand*(gbest(it,:)-particle(i).position);
%                             
%             particle(i).velocity=min(max(particle(i).velocity,-vmax),vmax);
%             
%             particle(i).position=particle(i).position+particle(i).velocity;
%             
%             particle(i).position=min(max(particle(i).position,xmin),xmax);
%             
%             particle(i).cost=Cost1(particle(i).position,Net);
%             
%             if particle(i).cost<particle(i).pbestcost
%                 particle(i).pbest=particle(i).position;
%                 particle(i).pbestcost=particle(i).cost;
% 
%                 if particle(i).pbestcost<gbestcost(it)
%                     gbest(it,:)=particle(i).pbest;
%                     gbestcost(it)=particle(i).pbestcost;
%                 end
%             end
%         end
%     end
%     
%     
%     w=w*wdamp;
%     
% end
% disp(['Iteration min Fc' num2str(oo) '/' num2str(numel(Fccy))]);
% TT(oo,:)=gbest(end,:);
% PSOResultsmin.Fc.result=TT;
% PSOResultsmin.Fc.iter{oo}=gbestcost;
% end
% %% Cu
% TT=zeros(numel(Cucy),nvar);
% 
% for oo=1:numel(Cucy)
% %%========================================
% w=1;
% wdamp=0.99;
% c1=2;
% c2=2;
% %========s0 mean  (kPa)=======
% xmin5=14.13;
% xmax5=400;
% dx5=xmax5-xmin5;
% vmax5=0.1*dx5;
% %======== Dr (%) =======
% xmin6=0; 
% xmax6=100;
% dx6=xmax6-xmin6;
% vmax6=0.1*dx6;
% %======== FC (%) =======
% xmin7=0; 
% xmax7=100;
% dx7=xmax7-xmin7;
% vmax7=0.1*dx7;
% %======== Cu =======
% xmin8=Cucy(oo);%1.52;
% xmax8=Cucy(oo);%28.12;
% dx8=xmax8-xmin8;
% vmax8=0.1*dx8;
% %======== D50 (mm) =======
% xmin9=0.029; 
% xmax9=0.46;
% dx9=xmax9-xmin9;
% vmax9=0.1*dx9;
% %============================cumulate all of restrictions ==============
% xmin=[xmin5;xmin6;xmin7;xmin8;xmin9];
% xmax=[xmax5;xmax6;xmax7;xmax8;xmax9];
% vmax=[vmax5;vmax6;vmax7;vmax8;vmax9];
% xmin=xmin';
% xmax=xmax';
% vmax=vmax';
% 
% empty_particle.position=[];
% empty_particle.velocity=[];
% empty_particle.cost=[];
% empty_particle.pbest=[];
% empty_particle.pbestcost=[];
% 
% particle=repmat(empty_particle,npop,1);
% 
% gbest=zeros(maxit,nvar);
% gbestcost=zeros(maxit,1);
% 
% for it=1:maxit
%     if it==1
%         gbestcost(1)=inf;
%         for i=1:npop
%             particle(i).velocity=zeros(1,nvar);
%             particle(i).position=xmin+(xmax-xmin).*rand(1,nvar);
%             particle(i).cost=Cost1(particle(i).position,Net);
%             particle(i).pbest=particle(i).position;
%             particle(i).pbestcost=particle(i).cost;
%             
%             if particle(i).pbestcost<gbestcost(it)
%                 gbest(it,:)=particle(i).pbest;
%                 gbestcost(it)=particle(i).pbestcost;
%             end
%         end
%     else
%         gbest(it,:)=gbest(it-1,:);
%         gbestcost(it)=gbestcost(it-1);
%         for i=1:npop
%             particle(i).velocity=w*particle(i).velocity...
%                                 +c1*rand*(particle(i).pbest-particle(i).position)...
%                                 +c2*rand*(gbest(it,:)-particle(i).position);
%                             
%             particle(i).velocity=min(max(particle(i).velocity,-vmax),vmax);
%             
%             particle(i).position=particle(i).position+particle(i).velocity;
%             
%             particle(i).position=min(max(particle(i).position,xmin),xmax);
%             
%             particle(i).cost=Cost1(particle(i).position,Net);
%             
%             if particle(i).cost<particle(i).pbestcost
%                 particle(i).pbest=particle(i).position;
%                 particle(i).pbestcost=particle(i).cost;
% 
%                 if particle(i).pbestcost<gbestcost(it)
%                     gbest(it,:)=particle(i).pbest;
%                     gbestcost(it)=particle(i).pbestcost;
%                 end
%             end
%         end
%     end
%       
%     w=w*wdamp;
%     
% end
% disp(['Iteration min Cu' num2str(oo) '/' num2str(numel(Cucy))]);
% TT(oo,:)=gbest(end,:);
% PSOResultsmin.Cu.result=TT;
% PSOResultsmin.Cu.iter{oo}=gbestcost;
% end
% %% D50
% TT=zeros(numel(D50cy),nvar);
% 
% for oo=1:numel(D50cy)
% %%========================================
% w=1;
% wdamp=0.99;
% c1=2;
% c2=2;
% %========s0 mean  (kPa)=======
% xmin5=14.13;
% xmax5=400;
% dx5=xmax5-xmin5;
% vmax5=0.1*dx5;
% %======== Dr (%) =======
% xmin6=0; 
% xmax6=100;
% dx6=xmax6-xmin6;
% vmax6=0.1*dx6;
% %======== FC (%) =======
% xmin7=0; 
% xmax7=100;
% dx7=xmax7-xmin7;
% vmax7=0.1*dx7;
% %======== Cu =======
% xmin8=1.52;
% xmax8=28.12;
% dx8=xmax8-xmin8;
% vmax8=0.1*dx8;
% %======== D50 (mm) =======
% xmin9=D50cy(oo);%0.029; 
% xmax9=D50cy(oo);%0.46;
% dx9=xmax9-xmin9;
% vmax9=0.1*dx9;
% %============================cumulate all of restrictions ==============
% xmin=[xmin5;xmin6;xmin7;xmin8;xmin9];
% xmax=[xmax5;xmax6;xmax7;xmax8;xmax9];
% vmax=[vmax5;vmax6;vmax7;vmax8;vmax9];
% xmin=xmin';
% xmax=xmax';
% vmax=vmax';
% 
% empty_particle.position=[];
% empty_particle.velocity=[];
% empty_particle.cost=[];
% empty_particle.pbest=[];
% empty_particle.pbestcost=[];
% 
% particle=repmat(empty_particle,npop,1);
% 
% gbest=zeros(maxit,nvar);
% gbestcost=zeros(maxit,1);
% 
% for it=1:maxit
%     if it==1
%         gbestcost(1)=inf;
%         for i=1:npop
%             particle(i).velocity=zeros(1,nvar);
%             particle(i).position=xmin+(xmax-xmin).*rand(1,nvar);
%             particle(i).cost=Cost1(particle(i).position,Net);
%             particle(i).pbest=particle(i).position;
%             particle(i).pbestcost=particle(i).cost;
%             
%             if particle(i).pbestcost<gbestcost(it)
%                 gbest(it,:)=particle(i).pbest;
%                 gbestcost(it)=particle(i).pbestcost;
%             end
%         end
%     else
%         gbest(it,:)=gbest(it-1,:);
%         gbestcost(it)=gbestcost(it-1);
%         for i=1:npop
%             particle(i).velocity=w*particle(i).velocity...
%                                 +c1*rand*(particle(i).pbest-particle(i).position)...
%                                 +c2*rand*(gbest(it,:)-particle(i).position);
%                             
%             particle(i).velocity=min(max(particle(i).velocity,-vmax),vmax);
%             
%             particle(i).position=particle(i).position+particle(i).velocity;
%             
%             particle(i).position=min(max(particle(i).position,xmin),xmax);
%             
%             particle(i).cost=Cost1(particle(i).position,Net);
%             
%             if particle(i).cost<particle(i).pbestcost
%                 particle(i).pbest=particle(i).position;
%                 particle(i).pbestcost=particle(i).cost;
% 
%                 if particle(i).pbestcost<gbestcost(it)
%                     gbest(it,:)=particle(i).pbest;
%                     gbestcost(it)=particle(i).pbestcost;
%                 end
%             end
%         end
%     end
%     
%     
%     w=w*wdamp;
%     
% end
% disp(['Iteration min D50' num2str(oo) '/' num2str(numel(D50cy))]);
% TT(oo,:)=gbest(end,:);
% PSOResultsmin.D50.result=TT;
% PSOResultsmin.D50.iter{oo}=gbestcost;
% end
%% All
%========s0 mean  (kPa)=======
xmin5=14.13;
xmax5=400;
dx5=xmax5-xmin5;
vmax5=0.1*dx5;
%======== Dr (%) =======
xmin6=0; 
xmax6=100;
dx6=xmax6-xmin6;
vmax6=0.1*dx6;
%======== FC (%) =======
xmin7=0; 
xmax7=100;
dx7=xmax7-xmin7;
vmax7=0.1*dx7;
%======== Cu =======
xmin8=1.52;
xmax8=28.12;
dx8=xmax8-xmin8;
vmax8=0.1*dx8;
%======== D50 (mm) =======
xmin9=0.029; 
xmax9=0.46;
dx9=xmax9-xmin9;
vmax9=0.1*dx9;
%============================cumulate all of restrictions ==============
xmin=[xmin5;xmin6;xmin7;xmin8;xmin9];
xmax=[xmax5;xmax6;xmax7;xmax8;xmax9];
vmax=[vmax5;vmax6;vmax7;vmax8;vmax9];
xmin=xmin';
xmax=xmax';
vmax=vmax';

empty_particle.position=[];
empty_particle.velocity=[];
empty_particle.cost=[];
empty_particle.pbest=[];
empty_particle.pbestcost=[];

particle=repmat(empty_particle,npop,1);

gbest=zeros(maxit,nvar);
gbestcost=zeros(maxit,1);

for it=1:maxit
    if it==1
        gbestcost(1)=inf;
        for i=1:npop
            particle(i).velocity=zeros(1,nvar);
            particle(i).position=xmin+(xmax-xmin).*rand(1,nvar);
            particle(i).cost=Cost1(particle(i).position,Net);
            particle(i).pbest=particle(i).position;
            particle(i).pbestcost=particle(i).cost;
            
            if particle(i).pbestcost<gbestcost(it)
                gbest(it,:)=particle(i).pbest;
                gbestcost(it)=particle(i).pbestcost;
            end
        end
    else
        gbest(it,:)=gbest(it-1,:);
        gbestcost(it)=gbestcost(it-1);
        for i=1:npop
            particle(i).velocity=w*particle(i).velocity...
                                +c1*rand*(particle(i).pbest-particle(i).position)...
                                +c2*rand*(gbest(it,:)-particle(i).position);
                            
            particle(i).velocity=min(max(particle(i).velocity,-vmax),vmax);
            
            particle(i).position=particle(i).position+particle(i).velocity;
            
            particle(i).position=min(max(particle(i).position,xmin),xmax);
            
            particle(i).cost=Cost1(particle(i).position,Net);
            
            if particle(i).cost<particle(i).pbestcost
                particle(i).pbest=particle(i).position;
                particle(i).pbestcost=particle(i).cost;

                if particle(i).pbestcost<gbestcost(it)
                    gbest(it,:)=particle(i).pbest;
                    gbestcost(it)=particle(i).pbestcost;
                end
            end
        end
    end
    
    
    w=w*wdamp;
    
end
TT=gbest(end,:);
PSOResultsmin.All.result=TT;
PSOResultsmin.All.iter=gbestcost;

%% maximization


%% S0
TT=zeros(numel(S0cy),nvar);

for oo=1:numel(S0cy)
%%========================================
w=1;
wdamp=0.99;
c1=2;
c2=2;
%========s0 mean  (kPa)=======
xmin5=S0cy(oo);%14.13;
xmax5=S0cy(oo);%400;
dx5=xmax5-xmin5;
vmax5=0.1*dx5;
%======== Dr (%) =======
xmin6=0; 
xmax6=100;
dx6=xmax6-xmin6;
vmax6=0.1*dx6;
%======== FC (%) =======
xmin7=0; 
xmax7=100;
dx7=xmax7-xmin7;
vmax7=0.1*dx7;
%======== Cu =======
xmin8=1.52;
xmax8=28.12;
dx8=xmax8-xmin8;
vmax8=0.1*dx8;
%======== D50 (mm) =======
xmin9=0.029; 
xmax9=0.46;
dx9=xmax9-xmin9;
vmax9=0.1*dx9;
%============================cumulate all of restrictions ==============
xmin=[xmin5;xmin6;xmin7;xmin8;xmin9];
xmax=[xmax5;xmax6;xmax7;xmax8;xmax9];
vmax=[vmax5;vmax6;vmax7;vmax8;vmax9];
xmin=xmin';
xmax=xmax';
vmax=vmax';

empty_particle.position=[];
empty_particle.velocity=[];
empty_particle.cost=[];
empty_particle.pbest=[];
empty_particle.pbestcost=[];

particle=repmat(empty_particle,npop,1);

gbest=zeros(maxit,nvar);
gbestcost=zeros(maxit,1);

for it=1:maxit
    if it==1
        gbestcost(1)=inf;
        for i=1:npop
            particle(i).velocity=zeros(1,nvar);
            particle(i).position=xmin+(xmax-xmin).*rand(1,nvar);
            particle(i).cost=Cost2(particle(i).position,Net);
            particle(i).pbest=particle(i).position;
            particle(i).pbestcost=particle(i).cost;
            
            if particle(i).pbestcost<gbestcost(it)
                gbest(it,:)=particle(i).pbest;
                gbestcost(it)=particle(i).pbestcost;
            end
        end
    else
        gbest(it,:)=gbest(it-1,:);
        gbestcost(it)=gbestcost(it-1);
        for i=1:npop
            particle(i).velocity=w*particle(i).velocity...
                                +c1*rand*(particle(i).pbest-particle(i).position)...
                                +c2*rand*(gbest(it,:)-particle(i).position);
                            
            particle(i).velocity=min(max(particle(i).velocity,-vmax),vmax);
            
            particle(i).position=particle(i).position+particle(i).velocity;
            
            particle(i).position=min(max(particle(i).position,xmin),xmax);
            
            particle(i).cost=Cost2(particle(i).position,Net);
            
            if particle(i).cost<particle(i).pbestcost
                particle(i).pbest=particle(i).position;
                particle(i).pbestcost=particle(i).cost;

                if particle(i).pbestcost<gbestcost(it)
                    gbest(it,:)=particle(i).pbest;
                    gbestcost(it)=particle(i).pbestcost;
                end
            end
        end
    end
    
    
    w=w*wdamp;
    
end
disp(['Iteration max S0' num2str(oo) '/' num2str(numel(S0cy))]);
TT(oo,:)=gbest(end,:);
PSOResultsmax.S0.result=TT;
PSOResultsmax.S0.iter{oo}=gbestcost;
end
% %% Dr
% TT=zeros(numel(Drcy),nvar);
% 
% for oo=1:numel(Drcy)
% %%========================================
% w=1;
% wdamp=0.99;
% c1=2;
% c2=2;
% %========s0 mean  (kPa)=======
% xmin5=14.13;
% xmax5=400;
% dx5=xmax5-xmin5;
% vmax5=0.1*dx5;
% %======== Dr (%) =======
% xmin6=Drcy(oo);%0; 
% xmax6=Drcy(oo);%100;
% dx6=xmax6-xmin6;
% vmax6=0.1*dx6;
% %======== FC (%) =======
% xmin7=0; 
% xmax7=100;
% dx7=xmax7-xmin7;
% vmax7=0.1*dx7;
% %======== Cu =======
% xmin8=1.52;
% xmax8=28.12;
% dx8=xmax8-xmin8;
% vmax8=0.1*dx8;
% %======== D50 (mm) =======
% xmin9=0.029; 
% xmax9=0.46;
% dx9=xmax9-xmin9;
% vmax9=0.1*dx9;
% %============================cumulate all of restrictions ==============
% xmin=[xmin5;xmin6;xmin7;xmin8;xmin9];
% xmax=[xmax5;xmax6;xmax7;xmax8;xmax9];
% vmax=[vmax5;vmax6;vmax7;vmax8;vmax9];
% xmin=xmin';
% xmax=xmax';
% vmax=vmax';
% 
% empty_particle.position=[];
% empty_particle.velocity=[];
% empty_particle.cost=[];
% empty_particle.pbest=[];
% empty_particle.pbestcost=[];
% 
% particle=repmat(empty_particle,npop,1);
% 
% gbest=zeros(maxit,nvar);
% gbestcost=zeros(maxit,1);
% 
% for it=1:maxit
%     if it==1
%         gbestcost(1)=inf;
%         for i=1:npop
%             particle(i).velocity=zeros(1,nvar);
%             particle(i).position=xmin+(xmax-xmin).*rand(1,nvar);
%             particle(i).cost=Cost2(particle(i).position,Net);
%             particle(i).pbest=particle(i).position;
%             particle(i).pbestcost=particle(i).cost;
%             
%             if particle(i).pbestcost<gbestcost(it)
%                 gbest(it,:)=particle(i).pbest;
%                 gbestcost(it)=particle(i).pbestcost;
%             end
%         end
%     else
%         gbest(it,:)=gbest(it-1,:);
%         gbestcost(it)=gbestcost(it-1);
%         for i=1:npop
%             particle(i).velocity=w*particle(i).velocity...
%                                 +c1*rand*(particle(i).pbest-particle(i).position)...
%                                 +c2*rand*(gbest(it,:)-particle(i).position);
%                             
%             particle(i).velocity=min(max(particle(i).velocity,-vmax),vmax);
%             
%             particle(i).position=particle(i).position+particle(i).velocity;
%             
%             particle(i).position=min(max(particle(i).position,xmin),xmax);
%             
%             particle(i).cost=Cost2(particle(i).position,Net);
%             
%             if particle(i).cost<particle(i).pbestcost
%                 particle(i).pbest=particle(i).position;
%                 particle(i).pbestcost=particle(i).cost;
% 
%                 if particle(i).pbestcost<gbestcost(it)
%                     gbest(it,:)=particle(i).pbest;
%                     gbestcost(it)=particle(i).pbestcost;
%                 end
%             end
%         end
%     end
%     
%     
%     w=w*wdamp;
%     
% end
% disp(['Iteration ' num2str(oo) '/' num2str(numel(Drcy))]);
% TT(oo,:)=gbest(end,:);
% PSOResultsmax.Dr.result=TT;
% PSOResultsmax.Dr.iter{oo}=gbestcost;
% end
% %% FC
% TT=zeros(numel(Fccy),nvar);
% 
% for oo=1:numel(Fccy)
% %%========================================
% w=1;
% wdamp=0.99;
% c1=2;
% c2=2;
% %========s0 mean  (kPa)=======
% xmin5=14.13;
% xmax5=400;
% dx5=xmax5-xmin5;
% vmax5=0.1*dx5;
% %======== Dr (%) =======
% xmin6=0; 
% xmax6=100;
% dx6=xmax6-xmin6;
% vmax6=0.1*dx6;
% %======== FC (%) =======
% xmin7=Fccy(oo);%0; 
% xmax7=Fccy(oo);%100;
% dx7=xmax7-xmin7;
% vmax7=0.1*dx7;
% %======== Cu =======
% xmin8=1.52;
% xmax8=28.12;
% dx8=xmax8-xmin8;
% vmax8=0.1*dx8;
% %======== D50 (mm) =======
% xmin9=0.029; 
% xmax9=0.46;
% dx9=xmax9-xmin9;
% vmax9=0.1*dx9;
% %============================cumulate all of restrictions ==============
% xmin=[xmin5;xmin6;xmin7;xmin8;xmin9];
% xmax=[xmax5;xmax6;xmax7;xmax8;xmax9];
% vmax=[vmax5;vmax6;vmax7;vmax8;vmax9];
% xmin=xmin';
% xmax=xmax';
% vmax=vmax';
% 
% empty_particle.position=[];
% empty_particle.velocity=[];
% empty_particle.cost=[];
% empty_particle.pbest=[];
% empty_particle.pbestcost=[];
% 
% particle=repmat(empty_particle,npop,1);
% 
% gbest=zeros(maxit,nvar);
% gbestcost=zeros(maxit,1);
% 
% for it=1:maxit
%     if it==1
%         gbestcost(1)=inf;
%         for i=1:npop
%             particle(i).velocity=zeros(1,nvar);
%             particle(i).position=xmin+(xmax-xmin).*rand(1,nvar);
%             particle(i).cost=Cost2(particle(i).position,Net);
%             particle(i).pbest=particle(i).position;
%             particle(i).pbestcost=particle(i).cost;
%             
%             if particle(i).pbestcost<gbestcost(it)
%                 gbest(it,:)=particle(i).pbest;
%                 gbestcost(it)=particle(i).pbestcost;
%             end
%         end
%     else
%         gbest(it,:)=gbest(it-1,:);
%         gbestcost(it)=gbestcost(it-1);
%         for i=1:npop
%             particle(i).velocity=w*particle(i).velocity...
%                                 +c1*rand*(particle(i).pbest-particle(i).position)...
%                                 +c2*rand*(gbest(it,:)-particle(i).position);
%                             
%             particle(i).velocity=min(max(particle(i).velocity,-vmax),vmax);
%             
%             particle(i).position=particle(i).position+particle(i).velocity;
%             
%             particle(i).position=min(max(particle(i).position,xmin),xmax);
%             
%             particle(i).cost=Cost2(particle(i).position,Net);
%             
%             if particle(i).cost<particle(i).pbestcost
%                 particle(i).pbest=particle(i).position;
%                 particle(i).pbestcost=particle(i).cost;
% 
%                 if particle(i).pbestcost<gbestcost(it)
%                     gbest(it,:)=particle(i).pbest;
%                     gbestcost(it)=particle(i).pbestcost;
%                 end
%             end
%         end
%     end
%     
%     
%     w=w*wdamp;
%     
% end
% disp(['Iteration max Fc' num2str(oo) '/' num2str(numel(Fccy))]);
% TT(oo,:)=gbest(end,:);
% PSOResultsmax.Fc.result=TT;
% PSOResultsmax.Fc.iter{oo}=gbestcost;
% end
% %% Cu
% TT=zeros(numel(Cucy),nvar);
% 
% for oo=1:numel(Cucy)
% %%========================================
% w=1;
% wdamp=0.99;
% c1=2;
% c2=2;
% %========s0 mean  (kPa)=======
% xmin5=14.13;
% xmax5=400;
% dx5=xmax5-xmin5;
% vmax5=0.1*dx5;
% %======== Dr (%) =======
% xmin6=0; 
% xmax6=100;
% dx6=xmax6-xmin6;
% vmax6=0.1*dx6;
% %======== FC (%) =======
% xmin7=0; 
% xmax7=100;
% dx7=xmax7-xmin7;
% vmax7=0.1*dx7;
% %======== Cu =======
% xmin8=Cucy(oo);%1.52;
% xmax8=Cucy(oo);%28.12;
% dx8=xmax8-xmin8;
% vmax8=0.1*dx8;
% %======== D50 (mm) =======
% xmin9=0.029; 
% xmax9=0.46;
% dx9=xmax9-xmin9;
% vmax9=0.1*dx9;
% %============================cumulate all of restrictions ==============
% xmin=[xmin5;xmin6;xmin7;xmin8;xmin9];
% xmax=[xmax5;xmax6;xmax7;xmax8;xmax9];
% vmax=[vmax5;vmax6;vmax7;vmax8;vmax9];
% xmin=xmin';
% xmax=xmax';
% vmax=vmax';
% 
% empty_particle.position=[];
% empty_particle.velocity=[];
% empty_particle.cost=[];
% empty_particle.pbest=[];
% empty_particle.pbestcost=[];
% 
% particle=repmat(empty_particle,npop,1);
% 
% gbest=zeros(maxit,nvar);
% gbestcost=zeros(maxit,1);
% 
% for it=1:maxit
%     if it==1
%         gbestcost(1)=inf;
%         for i=1:npop
%             particle(i).velocity=zeros(1,nvar);
%             particle(i).position=xmin+(xmax-xmin).*rand(1,nvar);
%             particle(i).cost=Cost2(particle(i).position,Net);
%             particle(i).pbest=particle(i).position;
%             particle(i).pbestcost=particle(i).cost;
%             
%             if particle(i).pbestcost<gbestcost(it)
%                 gbest(it,:)=particle(i).pbest;
%                 gbestcost(it)=particle(i).pbestcost;
%             end
%         end
%     else
%         gbest(it,:)=gbest(it-1,:);
%         gbestcost(it)=gbestcost(it-1);
%         for i=1:npop
%             particle(i).velocity=w*particle(i).velocity...
%                                 +c1*rand*(particle(i).pbest-particle(i).position)...
%                                 +c2*rand*(gbest(it,:)-particle(i).position);
%                             
%             particle(i).velocity=min(max(particle(i).velocity,-vmax),vmax);
%             
%             particle(i).position=particle(i).position+particle(i).velocity;
%             
%             particle(i).position=min(max(particle(i).position,xmin),xmax);
%             
%             particle(i).cost=Cost2(particle(i).position,Net);
%             
%             if particle(i).cost<particle(i).pbestcost
%                 particle(i).pbest=particle(i).position;
%                 particle(i).pbestcost=particle(i).cost;
% 
%                 if particle(i).pbestcost<gbestcost(it)
%                     gbest(it,:)=particle(i).pbest;
%                     gbestcost(it)=particle(i).pbestcost;
%                 end
%             end
%         end
%     end
%     
%   
%     w=w*wdamp;
%     
% end
% disp(['Iteration max Cu' num2str(oo) '/' num2str(numel(Cucy))]);
% TT(oo,:)=gbest(end,:);
% PSOResultsmax.Cu.result=TT;
% PSOResultsmax.Cu.iter{oo}=gbestcost;
% end
% %% D50
% TT=zeros(numel(D50cy),nvar);
% 
% for oo=1:numel(D50cy)
% %%========================================
% w=1;
% wdamp=0.99;
% c1=2;
% c2=2;
% %========s0 mean  (kPa)=======
% xmin5=14.13;
% xmax5=400;
% dx5=xmax5-xmin5;
% vmax5=0.1*dx5;
% %======== Dr (%) =======
% xmin6=0; 
% xmax6=100;
% dx6=xmax6-xmin6;
% vmax6=0.1*dx6;
% %======== FC (%) =======
% xmin7=0; 
% xmax7=100;
% dx7=xmax7-xmin7;
% vmax7=0.1*dx7;
% %======== Cu =======
% xmin8=1.52;
% xmax8=28.12;
% dx8=xmax8-xmin8;
% vmax8=0.1*dx8;
% %======== D50 (mm) =======
% xmin9=D50cy(oo);%0.029; 
% xmax9=D50cy(oo);%0.46;
% dx9=xmax9-xmin9;
% vmax9=0.1*dx9;
% %============================cumulate all of restrictions ==============
% xmin=[xmin5;xmin6;xmin7;xmin8;xmin9];
% xmax=[xmax5;xmax6;xmax7;xmax8;xmax9];
% vmax=[vmax5;vmax6;vmax7;vmax8;vmax9];
% xmin=xmin';
% xmax=xmax';
% vmax=vmax';
% 
% empty_particle.position=[];
% empty_particle.velocity=[];
% empty_particle.cost=[];
% empty_particle.pbest=[];
% empty_particle.pbestcost=[];
% 
% particle=repmat(empty_particle,npop,1);
% 
% gbest=zeros(maxit,nvar);
% gbestcost=zeros(maxit,1);
% 
% for it=1:maxit
%     if it==1
%         gbestcost(1)=inf;
%         for i=1:npop
%             particle(i).velocity=zeros(1,nvar);
%             particle(i).position=xmin+(xmax-xmin).*rand(1,nvar);
%             particle(i).cost=Cost2(particle(i).position,Net);
%             particle(i).pbest=particle(i).position;
%             particle(i).pbestcost=particle(i).cost;
%             
%             if particle(i).pbestcost<gbestcost(it)
%                 gbest(it,:)=particle(i).pbest;
%                 gbestcost(it)=particle(i).pbestcost;
%             end
%         end
%     else
%         gbest(it,:)=gbest(it-1,:);
%         gbestcost(it)=gbestcost(it-1);
%         for i=1:npop
%             particle(i).velocity=w*particle(i).velocity...
%                                 +c1*rand*(particle(i).pbest-particle(i).position)...
%                                 +c2*rand*(gbest(it,:)-particle(i).position);
%                             
%             particle(i).velocity=min(max(particle(i).velocity,-vmax),vmax);
%             
%             particle(i).position=particle(i).position+particle(i).velocity;
%             
%             particle(i).position=min(max(particle(i).position,xmin),xmax);
%             
%             particle(i).cost=Cost2(particle(i).position,Net);
%             
%             if particle(i).cost<particle(i).pbestcost
%                 particle(i).pbest=particle(i).position;
%                 particle(i).pbestcost=particle(i).cost;
% 
%                 if particle(i).pbestcost<gbestcost(it)
%                     gbest(it,:)=particle(i).pbest;
%                     gbestcost(it)=particle(i).pbestcost;
%                 end
%             end
%         end
%     end
%     
%     
%     
%     w=w*wdamp;
%     
% end
% disp(['Iteration max D50' num2str(oo) '/' num2str(numel(D50cy))]);
% TT(oo,:)=gbest(end,:);
% PSOResultsmax.D50.result=TT;
% PSOResultsmax.D50.iter{oo}=gbestcost;
% end
%% All
%========s0 mean  (kPa)=======
xmin5=14.13;
xmax5=400;
dx5=xmax5-xmin5;
vmax5=0.1*dx5;
%======== Dr (%) =======
xmin6=0; 
xmax6=100;
dx6=xmax6-xmin6;
vmax6=0.1*dx6;
%======== FC (%) =======
xmin7=0; 
xmax7=100;
dx7=xmax7-xmin7;
vmax7=0.1*dx7;
%======== Cu =======
xmin8=1.52;
xmax8=28.12;
dx8=xmax8-xmin8;
vmax8=0.1*dx8;
%======== D50 (mm) =======
xmin9=0.029; 
xmax9=0.46;
dx9=xmax9-xmin9;
vmax9=0.1*dx9;
%============================cumulate all of restrictions ==============
xmin=[xmin5;xmin6;xmin7;xmin8;xmin9];
xmax=[xmax5;xmax6;xmax7;xmax8;xmax9];
vmax=[vmax5;vmax6;vmax7;vmax8;vmax9];
xmin=xmin';
xmax=xmax';
vmax=vmax';

empty_particle.position=[];
empty_particle.velocity=[];
empty_particle.cost=[];
empty_particle.pbest=[];
empty_particle.pbestcost=[];

particle=repmat(empty_particle,npop,1);

gbest=zeros(maxit,nvar);
gbestcost=zeros(maxit,1);

for it=1:maxit
    if it==1
        gbestcost(1)=inf;
        for i=1:npop
            particle(i).velocity=zeros(1,nvar);
            particle(i).position=xmin+(xmax-xmin).*rand(1,nvar);
            particle(i).cost=Cost2(particle(i).position,Net);
            particle(i).pbest=particle(i).position;
            particle(i).pbestcost=particle(i).cost;
            
            if particle(i).pbestcost<gbestcost(it)
                gbest(it,:)=particle(i).pbest;
                gbestcost(it)=particle(i).pbestcost;
            end
        end
    else
        gbest(it,:)=gbest(it-1,:);
        gbestcost(it)=gbestcost(it-1);
        for i=1:npop
            particle(i).velocity=w*particle(i).velocity...
                                +c1*rand*(particle(i).pbest-particle(i).position)...
                                +c2*rand*(gbest(it,:)-particle(i).position);
                            
            particle(i).velocity=min(max(particle(i).velocity,-vmax),vmax);
            
            particle(i).position=particle(i).position+particle(i).velocity;
            
            particle(i).position=min(max(particle(i).position,xmin),xmax);
            
            particle(i).cost=Cost2(particle(i).position,Net);
            
            if particle(i).cost<particle(i).pbestcost
                particle(i).pbest=particle(i).position;
                particle(i).pbestcost=particle(i).cost;

                if particle(i).pbestcost<gbestcost(it)
                    gbest(it,:)=particle(i).pbest;
                    gbestcost(it)=particle(i).pbestcost;
                end
            end
        end
    end
    

    
    w=w*wdamp;
    
end
TT=gbest(end,:);
PSOResultsmax.All.result=TT;
PSOResultsmax.All.iter=gbestcost;
PSOResults{1,1}=PSOResultsmax;
PSOResults{1,2}=PSOResultsmin;