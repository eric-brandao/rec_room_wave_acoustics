function plot_acparameters_wa(general, geometry, sources, receivers, results, type)

if (strcmp(type,'decays')==1 || strcmp(type,'all')==1)
    %%%%% T20 %%%%%%%%%%%%%%%%%%%%%%%
    figure('Name','T20- Individuals')
    k=1;
    for js=1:length(sources)
        for jrec=1:length(receivers)
            T20=horzcat(results.T20.ind{js,jrec,:});
            parid=find(T20~=0); 
            semilogx(general.CenterFreq8va(parid),T20(parid),'o-','LineWidth',2); grid on; hold on;
            set(gca,'XTickLabel',{'16','31.5','63','125','250','500','1000','2000','4000','8000'},...
                'XTick',[16 31.5 63 125 250 500 1000 2000 4000 8000]);
            % set(gca,'Xticklabel',10.^get(gca,'Xtick'));
            xlabel('Frequency [Hz]'); ylabel('T20 [s]');
            legends{k}=strcat('source: ', num2str(js), ', receiver: ', num2str(jrec));
            k=k+1;
        end 
    end
    legends{k}='T20 mean and error';
       
    CI=[2*results.T20.sigma(parid);...
        2*results.T20.sigma(parid)];
    semilogx(general.CenterFreq8va(parid),results.T20.ux(parid),'o-k','LineWidth',2); grid on; hold on;
    
    if isfield(results,'decayMean')==1
        semilogx(general.CenterFreq8va(parid),results.T20.decayMean,...
            'd--','Color',[0 0.8 1],'LineWidth',2); grid on; hold on;
        k=k+1;
        legends{k}='T20 from mean decay'; hold on;        
    end
    shadedErrorBar(general.CenterFreq8va(parid), results.T20.ux(parid), CI'); hold on;
%     errorbar(general.CenterFreq8va,horzcat(results.T20.ux),2*results.T20.sigma,'o-k','LineWidth',2)
    k=k+1;
    legends{k}='Confidence interval';
    legend(legends,'Location','SouthWest');
    ylim([0 1.2*(max(results.T20.ux)+2*max(results.T20.sigma))]); xlim([16 1.2*general.CenterFreq8va(end)]);
    clear legends k
    
    %%%%% T30 %%%%%%%%%%%%%%%%%%%%%%%
    figure('Name','T30- Individuals')
    k=1;
    for js=1:length(sources)
        for jrec=1:length(receivers)
            T30=horzcat(results.T30.ind{js,jrec,:});
            semilogx(general.CenterFreq8va(parid),T30(parid),'o-','LineWidth',2); grid on; hold on;
            grid on; hold on;
            set(gca,'XTickLabel',{'16','31.5','63','125','250','500','1000','2000','4000','8000'},...
                'XTick',[16 31.5 63 125 250 500 1000 2000 4000 8000]);
            % set(gca,'Xticklabel',10.^get(gca,'Xtick'));
            xlabel('Frequency [Hz]'); ylabel('T30 [s]');
            legends{k}=strcat('source: ', num2str(js), ', receiver: ', num2str(jrec));
            k=k+1;
        end 
    end
    legends{k}='T30 mean and error';
    CI=[2*results.T30.sigma(parid);...
        2*results.T30.sigma(parid)];
    semilogx(general.CenterFreq8va(parid),results.T30.ux(parid),'o-k','LineWidth',2); grid on; hold on;
    if isfield(results,'decayMean')==1
        semilogx(general.CenterFreq8va(parid),results.T30.decayMean,...
            'd--','Color',[0 0.8 1],'LineWidth',2); grid on; hold on;
        k=k+1;
        legends{k}='T30 from mean decay'; hold on; 
    end
    shadedErrorBar(general.CenterFreq8va(parid), results.T30.ux(parid), CI'); hold on;    
%     errorbar(general.CenterFreq8va,horzcat(results.T30.ux),2*results.T30.sigma,'o-k','LineWidth',2)
    k=k+1;
    legends{k}='Confidence interval';
    legend(legends,'Location','SouthWest');
        ylim([0 1.2*(max(results.T30.ux)+2*max(results.T30.sigma))]); xlim([16 1.2*general.CenterFreq8va(end)]);
    
    clear legends k    
    %%%%% EDT %%%%%%%%%%%%%%%%%%%%%%%
    figure('Name','EDT - Individuals')
    k=1;
    for js=1:length(sources)
        for jrec=1:length(receivers)
            EDT=horzcat(results.EDT.ind{js,jrec,:});
            semilogx(general.CenterFreq8va(parid),EDT(parid),'o-','LineWidth',2); grid on; hold on;
            set(gca,'XTickLabel',{'16','31.5','63','125','250','500','1000','2000','4000','8000'},...
                'XTick',[16 31.5 63 125 250 500 1000 2000 4000 8000]);
            % set(gca,'Xticklabel',10.^get(gca,'Xtick'));
            xlabel('Frequency [Hz]'); ylabel('EDT [s]');
            legends{k}=strcat('source: ', num2str(js), ', receiver: ', num2str(jrec));
            k=k+1;
        end 
    end
    legends{k}='EDT mean and error';
    CI=[2*results.EDT.sigma(parid);...
        2*results.EDT.sigma(parid)];
    semilogx(general.CenterFreq8va(parid),results.EDT.ux(parid),'o-k','LineWidth',2); grid on; hold on;
    if isfield(results,'decayMean')==1
        semilogx(general.CenterFreq8va(parid),results.EDT.decayMean,...
            'd--','Color',[0 0.8 1],'LineWidth',2); grid on; hold on;
        k=k+1;
        legends{k}='EDT from mean decay';  hold on; 
    end
    shadedErrorBar(general.CenterFreq8va(parid), results.EDT.ux(parid), CI'); hold on;
%     errorbar(general.CenterFreq8va,horzcat(results.EDT.ux),2*results.EDT.sigma,'o-k','LineWidth',2)
    k=k+1;
    legends{k}='Confidence interval';
    legend(legends,'Location','SouthWest');
        ylim([0 1.2*(max(results.EDT.ux)+2*max(results.EDT.sigma))]); xlim([16 1.2*general.CenterFreq8va(end)]);
    clear legends k
end

%%
%%%%%%%%%%%%%%%%%%%%% C80 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (strcmp(type,'C80')==1 || strcmp(type,'c80')==1 || strcmp(type,'all')==1)
    figure('Name','C80')
    k=1;
    for js=1:length(sources)
        for jrec=1:length(receivers)
            C80=horzcat(results.C80.ind{js,jrec,:});
            semilogx(general.CenterFreq8va(parid),C80(parid),'o-','LineWidth',2); grid on; hold on;
            set(gca,'XTickLabel',{'16','31.5','63','125','250','500','1000','2000','4000','8000'},...
                'XTick',[16 31.5 63 125 250 500 1000 2000 4000 8000]);
            % set(gca,'Xticklabel',10.^get(gca,'Xtick'));
            xlabel('Frequency [Hz]'); ylabel('C_8_0 [s]');
            xlim([16 1.2*general.CenterFreq8va(end)]); 
            legends{k}=strcat('source: ', num2str(js), ', receiver: ', num2str(jrec));
            k=k+1;
        end
    end
    legends{k}='C80 mean and error';
    CI=[2*results.C80.sigma(parid);...
        2*results.C80.sigma(parid)];
    semilogx(general.CenterFreq8va(parid),results.C80.ux(parid),'o-k','LineWidth',2); grid on; hold on;
    shadedErrorBar(general.CenterFreq8va(parid), results.C80.ux(parid), CI'); hold on;
%     errorbar(general.CenterFreq8va,horzcat(results.C80.ux),2*results.C80.sigma,'o-k','LineWidth',2)
    k=k+1;
    legends{k}='Confidence interval';
    legend(legends,'Location','NorthWest');
    clear legends k
end
%%
%%%%%%%%%%%%%%%%%%%%% D50 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (strcmp(type,'D50')==1 || strcmp(type,'d50')==1 || strcmp(type,'all')==1)
    figure('Name','D50')
    k=1;
    for js=1:length(sources)
        for jrec=1:length(receivers)
            D50=horzcat(results.D50.ind{js,jrec,:});
            semilogx(general.CenterFreq8va(parid),D50(parid),'o-','LineWidth',2); grid on; hold on;
            set(gca,'XTickLabel',{'16','31.5','63','125','250','500','1000','2000','4000','8000'},...
                'XTick',[16 31.5 63 125 250 500 1000 2000 4000 8000]);
            % set(gca,'Xticklabel',10.^get(gca,'Xtick'));
            xlabel('Frequency [Hz]'); ylabel('D_5_0 [s]');
            ylim([0 1.0]); xlim([16 1.2*general.CenterFreq8va(end)]);
            legends{k}=strcat('source: ', num2str(js), ', receiver: ', num2str(jrec));
            k=k+1;
        end
    end
    legends{k}='D50 mean and error';
    CI=[2*results.D50.sigma(parid);...
        2*results.D50.sigma(parid)];
    semilogx(general.CenterFreq8va(parid),results.D50.ux(parid),'o-k','LineWidth',2); grid on; hold on;
    shadedErrorBar(general.CenterFreq8va(parid), results.D50.ux(parid), CI'); hold on;
    %     errorbar(general.CenterFreq8va,horzcat(results.D50.ux),2*results.D50.sigma,'o-k','LineWidth',2)
    k=k+1;
    legends{k}='Confidence interval';
    legend(legends,'Location','NorthWest');
    clear legends k
end

%%
%%%%%%%%%%%%%%%%%%%%% Ts %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (strcmp(type,'Ts')==1 || strcmp(type,'ts')==1 || strcmp(type,'all')==1)
    figure('Name','Ts')
    k=1;
    for js=1:length(sources)
        for jrec=1:length(receivers)
            Ts=horzcat(results.Ts.ind{js,jrec,:});
            semilogx(general.CenterFreq8va(parid),1000*Ts(parid),'o-','LineWidth',2); grid on; hold on;
            set(gca,'XTickLabel',{'16','31.5','63','125','250','500','1000','2000','4000','8000'},...
                'XTick',[16 31.5 63 125 250 500 1000 2000 4000 8000]);
            % set(gca,'Xticklabel',10.^get(gca,'Xtick'));
            xlabel('Frequency [Hz]'); ylabel('t_s [ms]');
            xlim([16 1.2*general.CenterFreq8va(end)]);
            legends{k}=strcat('source: ', num2str(js), ', receiver: ', num2str(jrec));
            k=k+1;
        end
    end
    legends{k}='Ts mean and error';
    CI=[2*results.Ts.sigma(parid);...
        2*results.Ts.sigma(parid)];
    semilogx(general.CenterFreq8va(parid),1000*results.Ts.ux(parid),'o-k','LineWidth',2); grid on; hold on;
    shadedErrorBar(general.CenterFreq8va(parid), 1000*results.Ts.ux(parid), 1000*CI'); hold on;
    %     errorbar(general.CenterFreq8va,1000*horzcat(results.Ts.ux),2000*results.Ts.sigma,'o-k','LineWidth',2)
    k=k+1;
    legends{k}='Confidence interval';
    legend(legends,'Location','SouthWest');
    clear legends k
end


% %%%%%%%%%%%%%%%%%%%%% G %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if (strcmp(type,'G')==1 || strcmp(type,'g')==1 || strcmp(type,'all')==1)
%     figure('Name','G')
%     k=1;
%     for js=sources.index
%         for jrec=receivers.pt.index
%             semilogx(general.CenterFreq8va,horzcat(results.raytracing.G.ind{js,jrec,:}),'o-','LineWidth',2); grid on; hold on;
%             set(gca,'XTickLabel',{'63','125','250','500','1000','2000','4000','8000'},...
%                 'XTick',[63 125 250 500 1000 2000 4000 8000]);
%             xlabel('Frequency [Hz]'); ylabel('G [dB]');
%             xlim([50 10000]); %ylim([0 1]);
%             legends{k}=strcat('source: ', num2str(js), ', receiver: ', num2str(jrec));
%             k=k+1;
%         end
%     end
%     legends{k}='G mean';
%     CI=[2*results.raytracing.G.sigma;...
%         2*results.raytracing.G.sigma];
%     semilogx(general.CenterFreq8va,horzcat(results.raytracing.G.ux),'o-k','LineWidth',2); grid on; hold on;
%     shadedErrorBar(general.CenterFreq8va, horzcat(results.raytracing.G.ux), CI'); hold on;
% %     semilogx(general.CenterFreq8va,horzcat(results.raytracing.G.ux),'o-k','LineWidth',2)
% %     errorbarlogx(0.03);
%     legend(legends,'Location','SouthWest');
% end