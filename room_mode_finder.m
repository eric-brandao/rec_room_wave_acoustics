%%%
% This function is used to find the modal frequencies of a rectangular room
% Inputs: Lx, Ly, Lz (room dimensions),
%         c0 - sound speed in m/s
%         fmax - target maximum frequency of analysis (the algorithm will
%         go beyond that for correct calculation of FRF
% Output: f_table (a table with each line containing the properties of a 
% modal frequency - collums mean:
%         1 - modal index
%         2 - modal frequency
%         3-5 - nx, ny, nz
%         6 - modal amplitude
%         7 - modal damping
%%%
function f_table=room_mode_finder(Lx,Ly,Lz,c0,fmax)

%% order
V=Lx*Ly*Lz;
od=2*ceil((2*fmax*V/c0)*sqrt(1./((Ly*Lz)^2+(Lx*Lz)^2+(Lx*Ly)^2)));
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(['Highest mode order is:' num2str(od)])
% od=ceil(2*fmax*Ly/c0);
% od=10;

%% indexes
nx=0:od; %length(nx)
ny=0:od; 
nz=0:od; 

% size_mesh=length(nx)*length(ny)*length(nz);
% hw = waitbar(0,'Encontrando modos...');
f_table=zeros(length(nx)*length(ny)*length(nz),5);
d=1;
for n=nz
    disp('%%%%%%%%%%%%% Be pacient... %%%%%%%%%%%%%%%%%')
    disp(['Number of modes to go is:' num2str(od^3-(n-1)*od^2)])
    disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

    for m=ny
        for l=nx
%             waitbar(d/size_mesh,hw)
            f(d)=(c0/2)*sqrt(((l/Lx)^2)+((m/Ly)^2)+((n/Lz)^2));
            f_table(d,:,:,:,:)=[d-1; f(d); l; m; n];
            d=d+1;
        end
    end
end
% close(hw)
f_table=f_table(2:end,:,:,:,:);

%% Associate mode energy and Tn*gamma
%%% For oblique modes
Axyz=(sqrt(1/8))*ones(length(f_table(:,1)),1);
enx=2; eny=2; enz=2;

Tn_gamma=(3*log(10)/(c0*((enx/Lx)+(eny/Ly)+(enz/Lz))))*ones(length(f_table(:,1)),1);


%%% For axial modes
axiidx=find(f_table(:,4)+f_table(:,5)==0); Axyz(axiidx)=sqrt(1/2);
enx=2; eny=1; enz=1;
Tn_gamma(axiidx)=(3*log(10)/(c0*((enx/Lx)+(eny/Ly)+(enz/Lz))));%*ones(length(Tn_gamma(axiidx)),1);%;

axiidy=find(f_table(:,3)+f_table(:,5)==0); Axyz(axiidy)=sqrt(1/2);
enx=1; eny=2; enz=1;
Tn_gamma(axiidy)=(3*log(10)/(c0*((enx/Lx)+(eny/Ly)+(enz/Lz))));%*ones(length(Tn_gamma(axiidy)),1);%;



axiidz=find(f_table(:,3)+f_table(:,4)==0); Axyz(axiidz)=sqrt(1/2);
enx=1; eny=1; enz=2;
Tn_gamma(axiidz)=(3*log(10)/(c0*((enx/Lx)+(eny/Ly)+(enz/Lz))));%*ones(length(Tn_gamma(axiidz)),1);%;



%%% For tangencial modes
tanidxy=find(f_table(:,3)~=0 & f_table(:,4)~=0 & f_table(:,5)==0); Axyz(tanidxy)=sqrt(1/4);
enx=2; eny=2; enz=1;
Tn_gamma(tanidxy)=(3*log(10)/(c0*((enx/Lx)+(eny/Ly)+(enz/Lz))));%*ones(length(Tn_gamma(tanidxy)),1);%;


tanidxz=find(f_table(:,3)~=0 & f_table(:,4)==0 & f_table(:,5)~=0); Axyz(tanidxz)=sqrt(1/4);
enx=2; eny=1; enz=2;
Tn_gamma(tanidxz)=(3*log(10)/(c0*((enx/Lx)+(eny/Ly)+(enz/Lz))));%*ones(length(Tn_gamma(tanidxz)),1);%;

tanidyz=find(f_table(:,3)==0 & f_table(:,4)~=0 & f_table(:,5)~=0); Axyz(tanidyz)=sqrt(1/4);
enx=1; eny=2; enz=2;
Tn_gamma(tanidyz)=(3*log(10)/(c0*((enx/Lx)+(eny/Ly)+(enz/Lz))));%*ones(length(Tn_gamma(tanidyz)),1);%;



% clc;
% Axyz
f_table=horzcat(f_table,Axyz,Tn_gamma);
disp(['Final frequency of table is: ' num2str(f_table(end,2)), ' [Hz]'])
disp(['There are ' num2str(f_table(end,1)), ' modes in the mode table'])
% f=sort(f);
% 
% [a,b]=hist(f,unique(f));