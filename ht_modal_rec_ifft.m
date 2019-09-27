%%%
% function to calculate the broad-band impulse response. Filters may work 
% only if you have the Audio Toolbox toolbox installed.
%%%

function [ht, decay, fc]=ht_modal_rec_ifft(general, geometry, sources, receivers,results)

%% Create time window if possible
if general.time(end) >= 1.2*geometry.T60
    len_ones=length(find(general.time<=1.2*geometry.T60));
    w_ones=ones(len_ones,1);
    len_hann=ceil((length(general.time)-len_ones)/4);
    w_hann=hanning(len_hann); w_hann=w_hann(len_hann/2:len_hann);
    w_zeros=zeros(length(general.time)-length(w_hann)-len_ones,1);
%     w_ones=ones(length(general.time)-length(w_hann),1);
    window=[w_ones; w_hann; w_zeros];
    clear len_ones len_hann w_hann w_ones w_zeros
else
    warning('The time response is shorter than 1.4 * T60. I can not apply a fade out to h(t)')
    window=ones(length(general.time),1);
end

%% Octave filter construction
try
    BW = '1 octave'; % Let's go on octave filtering
    N = 6;           % Filter Order
    fc = 1000;       % Center Frequency (Hz)
    % Fs = 48000;      % Sampling Frequency (Hz)
    oneOctaveFilter = octaveFilter('FilterOrder', N, ...
        'CenterFrequency', fc, 'Bandwidth', BW, 'SampleRate', general.Fs);
    fc = getANSICenterFrequencies(oneOctaveFilter);
    fc(fc<20) = [];
    fc(fc>general.fmax) = [];
    Nfc = length(fc);
    for j8va=1:Nfc
        fullOctaveFilterBank{j8va} = octaveFilter('FilterOrder', N, ...
            'CenterFrequency', fc(j8va), 'Bandwidth', BW, 'SampleRate', general.Fs); %#ok
    end

    %% Loop over all the sources and receivers
    ht(1:length(sources),1:length(receivers),1:Nfc+1)={zeros(1,length(general.time))};
    decay(1:length(sources),1:length(receivers),1:Nfc+1)={zeros(1,length(general.time))};
    for js=1:length(sources)
        for jrec=1:length(receivers)
            %%%% Step 1 - find the global impulse response and decay - it contains all
            %%%% frequencies - up to fmax
            ht(js,jrec,Nfc+1)={ifft(results.Hw{js,jrec},general.NFFT, 'symmetric')...
                .*window};
            decay(js,jrec,Nfc+1)={flipud(cumtrapz(general.time, flipud(abs(ht{js,jrec,Nfc+1}).^2)))}; %{cumsum(abs(results.ht{js,jrec}).^2,'reverse')};%
            %%%% Step 2 - Let's loop over the 1 octave frequency bands of analysis
            for j8va=1:Nfc
                oneOctaveFilter=fullOctaveFilterBank{j8va};
                ht(js,jrec,j8va)={oneOctaveFilter(ht{js,jrec,Nfc+1})};
                decay(js,jrec,j8va)={flipud(cumtrapz(general.time, flipud(abs(ht{js,jrec,j8va}).^2)))};
            end   
        end
    end
catch
    %% Loop over all the sources and receivers
    fc = 1;
    ht(1:length(sources),1:length(receivers))={zeros(1,length(general.time))};
    decay(1:length(sources),1:length(receivers))={zeros(1,length(general.time))};
    for js=1:length(sources)
        for jrec=1:length(receivers)
            %%%% Step 1 - find the global impulse response and decay - it contains all
            %%%% frequencies - up to fmax
            ht(js,jrec)={ifft(results.Hw{js,jrec},general.NFFT, 'symmetric')...
                .*window};
            decay(js,jrec)={flipud(cumtrapz(general.time, flipud(abs(ht{js,jrec}).^2)))}; %{cumsum(abs(results.ht{js,jrec}).^2,'reverse')};%
        end
    end
end %% end of try
