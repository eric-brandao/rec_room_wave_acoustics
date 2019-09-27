
function [ht, decay, decayMean]=ht_mode_recroom_meissner(general,geometry,sources,receivers,results)
% clc;
% Bn=geometry.damp; % dumping vs frequency for 20, 40, 80, 160, 320, 640, 1280, 2560;
% disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
% disp(['Damping is:' num2str(Bn)])
%% Room Volume
% V=Lx*Ly*Lz;
%% Calculate the modes;
% tic;
% f_table=room_mode_finder(geometry.Lx,geometry.Ly,geometry.Lz,general.c0,general.fmax);
% toc;
%% Calculete sound pressure at location x due to source at x0
% fn=results.fn_table(:,2);
% wn=2*pi*results.fn_table;
% kn=wn/general.c0;
% 
% nx=results.fn_table(:,3);
% ny=results.fn_table(:,4);
% nz=results.fn_table(:,5);
% Axyz=results.fn_table(:,6);
% 
% 
% %%% Meissner
% OmegaN=sqrt(wn.^2-(Bn)^2);

ht(1:length(sources),1:length(receivers),1:length(general.CenterFreq8va)+1)={[]};
decay(1:length(sources),1:length(receivers),1:length(general.CenterFreq8va)+1)={[]};
decayMean(1:length(general.CenterFreq8va))={zeros(1,length(general.time))};
h_all(1:length(sources),1:length(receivers))={zeros(1,length(general.time))};

for jf=1:length(general.CenterFreq8va)
    %%% Upper and lower frequencies in the 8va band
    flower = general.CenterFreq8va(jf) / (2^(0.5));
    fupper = general.CenterFreq8va(jf) * 2^(0.5);
    idf=find(results.fn_table(:,2)>=flower & results.fn_table(:,2)<fupper);
    
    fn=results.fn_table(idf,2);
    
    if isempty(fn)~=1
        nx=results.fn_table(idf,3);
        ny=results.fn_table(idf,4);
        nz=results.fn_table(idf,5);
        Axyz=results.fn_table(idf,6);
        gamma=(3*log(10)/geometry.T60)*(geometry.volume/(general.c0*geometry.St));
        T_n=results.fn_table(idf,7)/gamma;
        r_n=3*log(10)./T_n;
        
%         OmegaN=sqrt((2*pi*fn).^2-(Bn)^2);
        OmegaNN=sqrt((2*pi*fn).^2-(r_n).^2);
        
        %%% To calculate mean decay
        a_n=1./(2*r_n);
        b_n=r_n./(2*((r_n.^2)+(OmegaNN.^2)));
        c_n=OmegaNN./(2*((r_n.^2)+(OmegaNN.^2)));

        for js=1:length(sources)
            xs=sources(js).coord;
            Qs=sources(js).Q;

            %%% Variables independent of time
            psi_Xs=(Axyz.*(cos((nx*pi*xs(1))/geometry.Lx)).*(cos((ny*pi*xs(2))/geometry.Ly)).*...
                            (cos((nz*pi*xs(3))/geometry.Lz)))';


            for jrec=1:length(receivers)
            xr=receivers(jrec).coord;
            %%% Variables independent of time
            psi_Xr=Axyz.*(cos((nx*pi*xr(1))/geometry.Lx)).*(cos((ny*pi*xr(2))/geometry.Ly)).*...
                            (cos((nz*pi*xr(3))/geometry.Lz)); %% Collumn vector
            
            psi_M_xr=zeros(length(fn),length(general.time));
            
            if js==1 && jrec ==1
                psi_M_MeanDec=zeros(length(fn),length(general.time));
            end
            
            hw = waitbar(0,'1','Name','Where am I?');
            for jfn=1:length(fn)
                waitbar(jfn/length(fn),hw,strcat('Mounting matrix for source: ', num2str(js),...
                ', receiver: ', num2str(jrec),' and freq band: ', num2str(general.CenterFreq8va(jf))...
                ,' [Hz].'))
                %%% Assemble matrix
                %%% Constant damping over frequency
%                 psi_M_xr(jfn,:)=psi_Xr(jfn).*((exp(-Bn*general.time)).*...
%                     (sin(OmegaN(jfn)*general.time))./OmegaN(jfn)); %% Matrix fn x t
                %%% Damping varies with frequency
%                 psi_M_xr(jfn,:)=psi_Xr(jfn).*((exp(-r_n(jfn)*general.time)).*...
%                     (sin(OmegaNN(jfn)*general.time))./OmegaNN(jfn)); %% Matrix fn x t
                %%% Mount a matrix to use sum over modes
                psi_M_xr(jfn,:)=psi_Xr(jfn)*((exp(-r_n(jfn)*general.time)).*...
                    (sin(OmegaNN(jfn)*general.time))./OmegaNN(jfn))*psi_Xs(jfn); %% Matrix fn x t
                
                if js==1 && jrec ==1
                    psi_M_MeanDec(jfn,:)=((exp(-r_n(jfn)*general.time)).*...
                        (a_n(jfn)-b_n(jfn)*cos(2*OmegaNN(jfn)*general.time)+...
                        c_n(jfn)*sin(2*OmegaNN(jfn)*general.time)))/(OmegaNN(jfn)^2); %% Matrix fn x t
                end
            end
            delete(hw)
            %%% ht by matrix multiplication
%             ht(js,jrec,jf)={((general.rho0*Qs*general.c0^2)/geometry.volume)*...
%                 psi_Xs*psi_M_xr};
            %%% ht by sum of each collumn
            ht(js,jrec,jf)={((general.rho0*Qs*general.c0^2)/geometry.volume)*...
                sum(psi_M_xr,1)};
            decay(js,jrec,jf)={fliplr(cumtrapz(general.time, fliplr(abs(ht{js,jrec,jf}).^2)))};
            h_all(js,jrec)={h_all{js,jrec}+ht{js,jrec,jf}};
            ht(js,jrec,end)={h_all{js,jrec}};
            decay(js,jrec,end)={fliplr(cumtrapz(general.time, fliplr(abs(ht{js,jrec,end}).^2)))};
            end %% End of rec loop
        end %% End of source loop
    decayMean(jf)={(((general.rho0*general.c0^2)/geometry.volume)^2)*...
                sum(psi_M_MeanDec,1)};
    else %%%% In case there are no modes in that octave band.
        
        for js=1:length(sources)
            for jrec=1:length(receivers)
                ht(js,jrec)={[]};%{zeros(1,length(general.time))};
                decay(js,jrec)={[]};%{zeros(1,length(general.time))};
            end %% End of rec loop
        end %% End of source loop
    decayMean(jf)={[]};
    end %%% End if fn empty
    
end %% End of freq loop


% %%% Global response
% for js=1:length(sources)
%     for jrec=1:length(receivers)
%         h_all=vertcat(ht{js,jrec,:}); 
%         ht(jf,jrec,7)={sum(h_all)};
%     end %% End of rec loop
% end %% End of source loop