%%%
% function only used to plot the frequency response function for each 
% source - receiver pair
%%%
function plot_FRF(general, sources, receivers, results)

figure('Name','Magnitude of frequency response function','NumberTitle','off');
k=1;
for js=1:length(sources)
    for jrec=1:length(receivers)
        semilogx(general.freq,20*log10(abs(results.Hw{js,jrec})/sources(js).Q),'LineWidth',2); hold on;
        list{k}=strcat('source: ', num2str(js), '(Q = ', num2str(sources(js).Q), '[m^3/s]); receiver: ', num2str(jrec));
        k=k+1;
    end
end
grid on;
xlim([20 general.fmax]);  hold on; %ylim([-50 10]);
ax = gca;
set(ax,'XTick',[20 31 63 125 250 500 1000 2000 4000, 8000, 16000]); 
xlabel('Frequency [Hz]'); ylabel('Magnitude [dB]'); 
legend(list,'Location', 'NorthWest')

