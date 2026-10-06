function fh=PowerAnalysisVisualize(hz,avg_cort_spec,powerMap, global_sig_for, runInfo)
%Assumes HbT,Calcium, then calcium, or just

Cont2plot={'HbT'};
if ismember('Calcium',runInfo.Contrasts)
    Cont2plot={'HbT','Calcium'};
elseif ismember('FAD',runInfo.Contrasts)
    Cont2plot={'HbT','Calcium','FAD'};
end
            load(which('GoodWL.mat'),'xform_WL');
            colors={'r','b','k','g','c'};
            fh= figure('units','normalized','outerposition',[0 0 1 1]);          
            subplot(numel(Cont2plot)+1,2,1)            
            for contrast=1:numel(runInfo.Contrasts)
                if contrast<=3
                yyaxis left    
                loglog(hz,avg_cort_spec(:,contrast)./max(avg_cort_spec(:,contrast)),'Color',colors{contrast})               
                hold on
                ylabel('(\DeltaM)^2/Hz')
                else
                yyaxis right
                loglog(hz,avg_cort_spec(:,contrast)./max(avg_cort_spec(:,contrast)),'Color',colors{contrast})                       
                ylabel('(\DeltaF/F)^2/Hz')
                end
            end
            set(gca,'YScale','log')
            grid on;
            legend(runInfo.Contrasts,'Location','southwest')
            xlabel('Frequency (Hz)')
            title('Average Cortical Power Normalized by Max Power')

            subplot(numel(Cont2plot)+1,2,2)
                for contrast=1:numel(runInfo.Contrasts)
                if contrast<=3
                yyaxis left
                loglog(hz,global_sig_for(:,contrast)./max(global_sig_for(:,contrast)),'Color',colors{contrast})                       
                hold on
                ylabel('(\DeltaM)^2/Hz')
                else
                yyaxis right
                loglog(hz,global_sig_for(:,contrast)./max(global_sig_for(:,contrast)),'Color',colors{contrast})                    
                ylabel('(\DeltaF/F)^2/Hz')
                end
                end
            set(gca,'YScale','log')
            grid on;
            legend(runInfo.Contrasts,'Location','southwest')
            xlabel('Frequency (Hz)')
            title('Global Signal Power Normalized by Max Power')          

%Power Maps


ind=0;
for subplotz=3:2:2*numel(Cont2plot)+2
ind=ind+1;
            [~,contrastIndex]=ismember(Cont2plot{ind},runInfo.Contrasts);

            subplot(numel(Cont2plot)+1,2,subplotz)
            imagesc(powerMap(:,:,2,contrastIndex)) %isa
            axis image
            axis off
            c=colorbar;
            if contains(Cont2plot{1},'Hb')
            c.Label.String='(\DeltaM)^2/Hz'  ;
            else
            c.Label.String='(\DeltaF/F)^2/Hz'  ;                
            end
            colormap default
            caxis(quantile(reshape(powerMap(:,:,2,contrastIndex),128*128,[]),[0.05,.95]))  
            title(['ISA ', Cont2plot{ind}  ,' Power Map'])

            subplot(numel(Cont2plot)+1,2,subplotz+1)
            imagesc(powerMap(:,:,3,contrastIndex)) %delta
            axis off
            c=colorbar;
            axis image
            
            if contains(Cont2plot{1},'Hb')
            c.Label.String='(\DeltaM)^2/Hz'  ;
            else
            c.Label.String='(\DeltaF/F)^2/Hz'   ;               
            end          
            colormap default
            caxis(quantile(reshape(powerMap(:,:,3,contrastIndex),128*128,[]),[0.05,.95]))  
            title(['Delta ', Cont2plot{ind}  ,' Power Map'])            
    
end


end