%%% Author: Eric Brandão
%%% Last significant modification: 26/12/2018

%% Description
%%% This function calculates the acoustical objective parameters obtained
%%% from impulse responses
%%% Also, the travelling time from source to receiver has been eliminated
%%% from relevant calculations, so that exactness may be improved.

%% Input variables
%%%% general - Geometrical information of the room - struct with planes and normals
%%%% sources - to have the number of receivers
%%%% receivers - to have the number of receivers
%%%% results
% %%%% time - time of arrival of energies for each set of source-receiver-frequency
% %%%% energy - energy of arrival of energies for each set of source-receiver-frequency
% %%%% decay - decay curve for each set of source-receiver-frequency


function [T20,T30,EDT,C80,D50,Ts]=ac_parameters_ht(general,sources,receivers,results)


%%% Find out at which index the direct sound is
% if strcmp(method,'reflecto')==1
%     [~,tDir_id] = findpeaks(results.raytracing.reflecto{1,1,1});
%     tDir_id=tDir_id(1);
% end


for js=1:length(sources) %% loop over sources
    for jrec=1:length(receivers) %% loop over receivers
        for jf=1:length(general.CenterFreq8va)  %% Loop over frequencies

            %%%% Estimate time of arrival of the direct sound in the
            %%%% reflectogram context (it may be an histogram)
            tDir=0.9*norm(sources(js).coord-receivers(jrec).coord)/general.c0;
            tDir_id=find(general.time<=tDir); tDir_id=tDir_id(end);
            
%             [~,tDir_id] = findpeaks(results.ht{js,jrec,1});
%             tDir_id=tDir_id(1);
                        
            ht=results.ht{js,jrec,jf}; %% the reflectogram
            
            if isempty(ht)~=1
                if isrow(ht)~=1
                    ht=(ht(tDir_id:end))';
                    dec=results.decay{js,jrec,jf};
                    dec=(dec(tDir_id:end))';
                else
                    ht=ht(tDir_id:end);
                    dec=results.decay{js,jrec,jf}; 
                    dec=dec(tDir_id:end);
                end
                
                
                t=general.time;
                t=t(tDir_id:end); 
                %% T20, T30, EDT

                decaydB=10*log10(dec./max(dec));
                a10=find(decaydB<-0.005 & decaydB > -10.005); 
                a20=find(decaydB<-5 & decaydB > -5-20); 
                a30=find(decaydB<-5 & decaydB > -5-30);
                %%% EDT
                idi10=a10(1); idf10=a10(end);
                p10 = polyfit(t(idi10:idf10),decaydB(idi10:idf10),1);
                htFIT10=p10(1)*t+p10(2);
                EDT.ind(js,jrec,jf)={-60/p10(1)};

                %%% T20, T30
                idi20=a20(1); idf20=a20(end);
                idi30=a30(1); idf30=a30(end);
                p20 = polyfit(t(idi20:idf20),decaydB(idi20:idf20),1); 
                p30 = polyfit(t(idi30:idf30),decaydB(idi30:idf30),1);

                htFIT20=p20(1)*t+p20(2);
                htFIT30=p30(1)*t+p30(2);

                T20.ind(js,jrec,jf)={-60/p20(1)};
                T30.ind(js,jrec,jf)={-60/p30(1)};
                
                %% C80, D50, Ts
                id80=find(t<=80/1000+t(1)); 
                C80.ind(js,jrec,jf)={10*log10(trapz(t(id80),ht(id80).^2)/...
                    trapz(t(id80(end)+1:end),ht(id80(end)+1:end).^2))};

                id50=find(t<=50/1000+t(1)); 
                D50.ind(js,jrec,jf)={trapz(t(id50),ht(id50).^2)/trapz(t(1:end),ht(1:end).^2)};

                Ts.ind(js,jrec,jf)={trapz(t,t.*ht.^2)/trapz(t,ht.^2)};
    %             %% G
    % %             G.ind(js,jrec,jf)={10*log10(sum(ref))-10*log10(sources.Wlin(js,jf)/(4*pi*10^2))};
    %             G.ind(js,jrec,jf)={10*log10(trapz(ref))-10*log10(sources.Wlin(js,jf)/(4*pi*10^2))};
            else
                T20.ind(js,jrec,jf)={0};
                T30.ind(js,jrec,jf)={0};
                EDT.ind(js,jrec,jf)={0};
                C80.ind(js,jrec,jf)={0};
                D50.ind(js,jrec,jf)={0};
                Ts.ind(js,jrec,jf)={0};
            end
        end
    end
end
% Mean ans STD values
d=1;

       
for js=1:length(sources) %% loop over sources
    for jrec=1:length(receivers) %% loop over receivers
        T20list(d,:)=horzcat(T20.ind{js,jrec,:});
        T30list(d,:)=horzcat(T30.ind{js,jrec,:});
        EDTlist(d,:)=horzcat(EDT.ind{js,jrec,:});
        C80list(d,:)=horzcat(C80.ind{js,jrec,:});
        D50list(d,:)=horzcat(D50.ind{js,jrec,:});
        Tslist(d,:)=horzcat(Ts.ind{js,jrec,:});
%         Glist(d,:)=horzcat(G.ind{js,jrec,:});
        d=d+1;
    end
end

T20.ux=mean(T20list,1);
T30.ux=mean(T30list,1);
EDT.ux=mean(EDTlist,1);
C80.ux=mean(C80list,1);
D50.ux=mean(D50list,1);
Ts.ux=mean(Tslist,1);
% G.ux=mean(Glist,1);

T20.sigma=std(T20list,1,1);
T30.sigma=std(T30list,1,1);
EDT.sigma=std(EDTlist,1,1);
C80.sigma=std(C80list,1,1);
D50.sigma=std(D50list,1,1);
Ts.sigma=std(Tslist,1,1);
% G.sigma=std(Glist,1);
%% From mean decays
if isfield(results,'decayMean')==1 
    d=1;
    for jf=1:length(general.CenterFreq8va)  %% Loop over frequencies
        if isempty(results.decayMean{jf}) ~=1
            t=general.time;
            decM=results.decay{1,1,jf};
            decaydB=10*log10(decM./max(decM));
            a10=find(decaydB<-0.005 & decaydB > -10.005); 
            a20=find(decaydB<-5 & decaydB > -5-20); 
            a30=find(decaydB<-5 & decaydB > -5-30);
            idi10=a10(1); idf10=a10(end);
            p10 = polyfit(t(idi10:idf10),decaydB(idi10:idf10),1);
            htFIT10=p10(1)*t+p10(2);
            EDT.decayMean(d)=(-60-p10(2))/p10(1);
            
            %%% T20, T30
            idi20=a20(1); idf20=a20(end);
            idi30=a30(1); idf30=a30(end);
            p20 = polyfit(t(idi20:idf20),decaydB(idi20:idf20),1); 
            p30 = polyfit(t(idi30:idf30),decaydB(idi30:idf30),1);

            htFIT20=p20(1)*t+p20(2);
            htFIT30=p30(1)*t+p30(2);

            T20.decayMean(d)=(-60-p20(2))/p20(1);
            T30.decayMean(d)=(-60-p30(2))/p30(1);
            
            d=d+1;
        end
    end 
end

                %%%% Of mean decays in case it existis
%                 if isfield(results,'decayMean')==1 && js==1 && jrec==1
%                     if isempty(results.decayMean{jf}) ~= 1
%                     end
%                 end



% % plot(time,htFIT10); hold on;
% % plot(time(idi10:idf10),decaydB(idi10:idf10))
% 
% 
