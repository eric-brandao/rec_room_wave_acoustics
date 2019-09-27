%%%
% Function to plot impulse response and decay
%%%
function plot_ht_decay(general, sources, receivers, results, colorplt, range)
if general.CenterFreq8va == 1
    Nfc = 1;
else
    Nfc=length(general.CenterFreq8va)+1;
end
%% colors
if colorplt==0
    ColorVec1=[0 0 1]; %% blue
    ColorVec2=[1 0 0.6]; %% magenta
elseif colorplt==1
    ColorVec1=[0.5 0.5 0.5]; %% grey
    ColorVec2=[0 0 0]; %% black
    
else
    ColorVec1=[1 0 0]; %% red
    ColorVec2=[0.2 0.6 0.4]; %% green
end
%% range
if range < 10
    range=-range;
elseif range==10
    range=80;
end

if range+10<20
    range=10;
    warning('The specified range is too small. I am plotting 10 [dB] of decay.')
end
%% Impulse response and decay
for js=1:length(sources)
    for jrec=1:length(receivers)    
        figure('Name',strcat('Impulse response and decay (dB): source: ', num2str(js), '; and receiver: ',...
            num2str(jrec), '.'),'units','normalized','outerposition',[0 0 1 1])
        d=1;
        for j8va=1:Nfc
            if Nfc == 1
                subplot(1,1,d)
            elseif Nfc>1 && Nfc<=3
                subplot(1,3,d)
            elseif Nfc>3 && Nfc <=8
                subplot(2,4,d)
            else
                subplot(3,4,d)
            end
            if isempty(results.ht{js,jrec,j8va})~=1
                plot(general.time,10*log10(abs(results.ht{js,jrec,j8va}/...
                    max(results.ht{js,jrec,j8va})).^2),'LineWidth',1,'Color',ColorVec1); hold on;
                plot(general.time,10*log10(results.decay{js,jrec,j8va}/...
                    max(results.decay{js,jrec,j8va})),'LineWidth',2,'Color',ColorVec2); hold on;
%                 if isfield(results,'decayMean')==1 && j8va ~= Nfc(end)
%                     plot(general.time,10*log10(abs(results.decayMean{j8va})/...
%                         max(results.decayMean{j8va})),'LineWidth',2); hold on;
%                 end
            end
            d=d+1;
            ylim([-range 10])
            if j8va ~= Nfc(end)
                title(strcat(num2str(general.CenterFreq8va(j8va)), ' [Hz]'))
            else
                title('Global')
            end
        end %% End 8va bands
    end %% End receivers
end %% End sources


%% Impulse response
for js=1:length(sources)
    for jrec=1:length(receivers)    
        figure('Name',strcat('Impulse response (linear): source: ', num2str(js), '; and receiver: ',...
            num2str(jrec), '.'),'units','normalized','outerposition',[0 0 1 1])
        d=1;
        for j8va=1:Nfc
            if Nfc == 1
                subplot(1,1,d)
            elseif Nfc>1 && Nfc<=3
                subplot(1,3,d)
            elseif Nfc>3 && Nfc <=8
                subplot(2,4,d)
            else
                subplot(3,4,d)
            end
            if isempty(results.ht{js,jrec,j8va})~=1
                plot(general.time,results.ht{js,jrec,j8va}/...
                    max(results.ht{js,jrec,j8va}),'b','LineWidth',1,'Color',ColorVec1); hold on;
            end
            d=d+1;
            ylim([-1.1 1.1])
            if j8va ~= Nfc(end)
                title(strcat(num2str(general.CenterFreq8va(j8va)), ' [Hz]'))
            else
                title('Global')
            end
        end %% End 8va bands
    end %% End receivers
end %% End sources
