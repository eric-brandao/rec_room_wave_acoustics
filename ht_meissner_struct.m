clc; clear all; close all
%%% Para calcular a resposta de um sistema fonte-receptor e uma sala retangular de paredes 
%%% rígidas (com amortecimento). O cálculo é feito no domínio do tempo, de acordo com Meissner.  
%%% Neste caso, para a resposta ao impulso (domínio do tempo) é
%%% obtida diretamente por banda de oitava não sendo necessário filtra-la.
%%% Para calcular o espectro aplica-se a FFT na h(t).
%%% As curvas de decaimento podem ser obtidas de duas maneiras: (1) Pelo
%%% método descrito por Meissner (o que implica em um considerável custo
%%% computacional); (2) Pela integração reversa de Schroeder.
%% General - Air properties
general.Temp=20; %% temperature [ºC]
general.HR=50; %% relative umidity [%]
general.P0=101325; %% Atmospheric pressure [Pa]
[general.rho0,general.c0]=propair_panneton(general.Temp,general.P0,general.HR);
% general.rho0=1.21;
% general.c0=343;
%% General Algorithm Properties
%%%% User data %%%%%%%%%%%%%%%
general.fmax=500; %% Maximum target octave band
general.CenterFreq8va=(10^3)*(2.^(-6:4)); 
general.CenterFreq8va=general.CenterFreq8va(find(general.CenterFreq8va<=general.fmax));
general.order=13; %% Order of the FFT
%%%%%% Calculated %%%%%%%%%%%%%
general.Fs=6000; %% Sample rate > 2*highest frequency to run
general.Dt=1/general.Fs; %% time discretization
general.NFFT=2^(general.order); %% Number of points in the FFT

%%% The time vector can be created in the following way. It is more similar
%%% to frequency dommain calculations. But, it can also be created more
%%% freely (i.e. time=0:1/Fs:1.2). In that case, the number of points on
%%% the frequency vector (and all spectrum) will have to be done
%%% accordingly.
general.time=linspace(0,(general.NFFT-1)*general.Dt,general.NFFT);
general.freq=linspace(0,general.Fs/2,general.NFFT/2); %% Frequency vector
general.Df=general.freq(2)-general.freq(1); %% Frequency discretization

% general.idf=find(general.freq<=1.1*general.fmax);
% general.idf=general.idf(end)+1;
disp('******** Signal parameters ********')
disp(['Frequency Resolution=' num2str(general.Df) ' [Hz]'])
disp(['Total time of impulse response=' num2str(general.time(end)) ' [s]'])

%% Room geometry
%%%% User data %%%%%%%%%%%%%%
geometry.Lx=1.54*4;   % 8;%  %%% Length 4.06;%
geometry.Ly=1.28*4;   % 5;%  %%% Width 3.36;%
geometry.Lz=4; %%%%     3;%      Heigth 2.8;%
geometry.T60=1;%2; %%% Reverberation time (Global or per-frequency)

%%%% Calculated data %%%%%%%%%%%%%%%
geometry.volume=geometry.Lx*geometry.Ly*geometry.Lz; %%% Room volume
geometry.freqSchroeder=2000*sqrt(max(geometry.T60)/geometry.volume); %% Schroeder frequency based on max T60 (conservative estimative)
geometry.St=2*geometry.Lx*geometry.Ly+2*geometry.Lx*geometry.Lz+2*geometry.Ly*geometry.Lz;
% geometry.alpha=1-exp((geometry.volume/(geometry.St*geometry.T60))*(4*log(10e-6)/general.c0));
geometry.damp=3*log(10)/geometry.T60;

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(['Schroeder Frequency is:' num2str(geometry.freqSchroeder)])
%% Sources
sources(1).coord=[0.01 0.01 0.01]; %[2 3 1];%%% Source coordinates
sources(1).Q=1;%1; %%% Source Volume velocity (may be a scalar or a vector same size as freq)
% sources(2).coord=[1 1 1]; 
% sources(2).Q=0.8;
% sources.Ns=length(sources);

%% Receivers
receivers(1).coord=[geometry.Lx-0.01 geometry.Ly-0.01 geometry.Lz-0.01];%[7 4 1.5];%[4 2 1.5];%
receivers(2).coord=[3 3 1.2];
receivers(3).coord=[geometry.Lx/2 geometry.Ly/2 1.2];

%% Fill the table of acoustic modes
results.fn_table=room_mode_finder(geometry.Lx,geometry.Ly,geometry.Lz,general.c0,general.fmax);

%% Fill results struct with empty data to pass to the mex function
% for js=1:length(sources)
%     for jrec=1:length(receivers)
%         results.Hw(js,jrec)={[]};
%     end
% end
% clear js jrec
%% Call Mex to calculate h(t)
tic;
[results.ht, results.decay, results.decayMean]=ht_mode_recroom_meissner(general,geometry,sources,receivers,results);
toc;
%% Calculate FRF (global and octave filtered)
for js=1:length(sources) %% loop over sources
    for jrec=1:length(receivers) %% loop over receivers
        Hww=2*fft(results.ht{js,jrec,end})/(general.NFFT); Hww=Hww(1:general.NFFT/2);
        results.Hw(js,jrec)={2*pi*general.freq.*Hww};
    end
end
clear Hww js jrec
%% Calculate acoustical parameters
[results.T20,results.T30,results.EDT,results.C80,results.D50,results.Ts]=...
    ac_parameters_ht(general,sources,receivers,results);

%% plots
plot_FRF(general, sources, receivers, results)
plot_ht_decay(general, sources, receivers, results,0,35)
plot_acparameters_wa(general, geometry, sources, receivers, results, 'all')




