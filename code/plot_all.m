%% Make Plots

% Plot the raw evoked activity correlation
plot_raw_corr(corrMatrix)

% Plot top noise components by stimulus
plot_ncs_evecs(evecs, evals)

% Plot the evoked activity, zscored over time, by neuron and condition
plot_zscored_evoked(evoked,stim_trials)

% Plot 3 dim noise correlation graphs
plot_3d_ncs_on_sig(sig, ncs)

% Plot the noise correlation matrices
plot_nc_mats(ncs, ranks1, ranks2, ranks3)

% Plot the noise corrs as functions of signal tuning
plot_sig_vs_noise(sigcorr, ncs, ncs_avg)

% Plot the signal correlation matrix as points in signal space
plot_normed_sigcorr(sig_normed)

% Plot the rankings of how strongly each neurons correlates with others by signal
plot_total_corr_ranks(ranks1,ranks2,ranks3)


%% Plotting functions

function [] = plot_raw_corr(corrMatrix)

   figure()
   imagesc(corrMatrix)
   xlabel('Cell Number')
   ylabel('Cell Number')
   title('Raw Correlation Across Cells')
   colorbar()
end

function [] = plot_ncs_evecs(evecs, evals)

   figure()
   set(gcf,'Position', [95, 528, 1520, 416])
   colorOrder = get(gca, 'ColorOrder');
   subplot(1,2,1)
   plot(evecs{1}(:,1), '-o', 'Color', colorOrder(1,:)); hold on
   plot(evecs{2}(:,1), '-o', 'Color', colorOrder(2,:));
   plot(evecs{3}(:,1), '-o', 'Color', colorOrder(3,:));

   plot( evecs{1}(:,2), '--o', 'Color', colorOrder(1,:));
   plot( evecs{2}(:,2), '--o', 'Color', colorOrder(2,:));
   plot(-evecs{3}(:,2), '--o', 'Color', colorOrder(3,:));

   plot([1,11], [0,0],'--k')
   xlim([1,11])
   grid on

   title('Noise Correlation Matrix Eigenvectors')
   xlabel('Neuron ID')
   ylabel('Loading')
   legend({'S1 E1', 'S2 E1', 'S3 E1', 'S1 E2', 'S2 E2', 'S3 E2'},'Location','NorthWest')

   subplot(1,2,2)
   plot(evals{1}, '-o', 'Color', colorOrder(1,:)); hold on
   plot(evals{2}, '-o', 'Color', colorOrder(2,:));
   plot(evals{3}, '-o', 'Color', colorOrder(3,:));
   grid on

   title('Noise Correlation Matrix Spectra')
   xlabel('Eigenvalue Rank')
   ylabel('Eigenvalue')
   legend({'S1', 'S2', 'S3'})
   xlim([1,11])
end

function [] = plot_zscored_evoked(evoked, stim_trials)
   figure()
   
   for i = 1:3
      subplot(3,1,i)
      
      z_evoked = zscore(evoked(:,stim_trials{i})');
      z_evoked_mu = mean(z_evoked,2);
      z_evoked_sd = std(z_evoked,0,2);

      ax = plot(z_evoked_mu); hold on
      plot(z_evoked_mu + z_evoked_sd, '--', 'Color', ax.Color)
      plot(z_evoked_mu - z_evoked_sd, '--', 'Color', ax.Color)
      grid on
      
      title(['S', num2str(i)])
      xlabel('Trial Number')
      ylabel('Z-Scored dF/F')
      legend({'mean', 'sd'})
      xlim([1,length(z_evoked_mu)])
      ylim([-3,3])
   end
end

function [] = plot_3d_ncs_on_sig(sig, ncs)
   figure()
   set(gcf,'Position', [59         491        1533         374])
   for k = 1:3
      subplot(1,3,k)

      scatter3(sig(1,:), sig(2,:), sig(3,:))
      hold on
      zlim([0,60])
      ylim([0,60])
      xlim([0,60])

      title(['Noise corr. graph given S', num2str(k)])
      xlabel('S1 response')
      ylabel('S2 resposne')
      zlabel('S3 response')


      minNC = min(   get_upper(ncs{k}));
      medNC = median(get_upper(ncs{k}));
      maxNC = max(   get_upper(ncs{k}));

      lcolors = jet(101);
      for i = 1:10
         for j = (i+1):11
            norm_val(i,j) = (ncs{k}(i,j)-minNC)/(maxNC-minNC);
            link_col{i,j} = lcolors(round(norm_val(i,j)*100)+1,:);
         end
      end

      quad4 = prctile(get_upper(norm_val*100),0);
      for i = 1:10
         for j = 2:11
            if norm_val(i,j)*100 >= quad4
               plot3(sig(1,[i,j]), sig(2,[i,j]), sig(3,[i,j]), 'Color', link_col{i,j}, 'LineWidth', 2);
            end
         end
      end

   end
end


function [] = plot_nc_mats(ncs, ranks1, ranks2, ranks3)

   figure()

   subplot(3,3,1)
   imagesc(ncs{1})
   xlabel('Neuron ID')
   ylabel('Neuron ID')
   title('S1 NCs')

   subplot(3,3,2)
   imagesc(ncs{2})
   xlabel('Neuron ID')
   ylabel('Neuron ID')
   title('S2 NCs')

   subplot(3,3,3)
   imagesc(ncs{3})
   xlabel('Neuron ID')
   ylabel('Neuron ID')
   title('S3 NCs')

   subplot(3,3,4)
   imagesc(ncs{1}(ranks1,ranks1))
   xlabel('Neuron Total Corr. Rank')
   ylabel('Neuron Total Corr. Rank')
   title('S1 NCs')

   subplot(3,3,5)
   imagesc(ncs{2}(ranks2,ranks2))
   xlabel('Neuron Total Corr. Rank')
   ylabel('Neuron Total Corr. Rank')
   title('S2 NCs')

   subplot(3,3,6)
   imagesc(ncs{3}(ranks3,ranks3))
   xlabel('Neuron Total Corr. Rank')
   ylabel('Neuron Total Corr. Rank')
   title('S3 NCs')

   subplot(3,3,7)
   histogram(get_upper(ncs{1}),'BinEdges',-1:0.05:1)
   grid on
   title('S1 NC Histogram')
   ylabel('Count')
   xlabel('Noise Correlation')

   subplot(3,3,8)
   histogram(get_upper(ncs{2}),'BinEdges',-1:0.05:1)
   grid on
   title('S2 NC Histogram')
   ylabel('Count')
   xlabel('Noise Correlation')

   subplot(3,3,9)
   histogram(get_upper(ncs{3}),'BinEdges',-1:0.05:1)
   grid on
   title('S3 NC Histogram')
   ylabel('Count')
   xlabel('Noise Correlation')

end

function [] = plot_sig_vs_noise(sigcorr, ncs, ncs_avg)

   sigvec = get_upper(sigcorr);

   figure()
   for i = 1:3
      subplot(2,2,i)
      
      ncsvec = get_upper(ncs{i});
      [m, b] = fit_orthog_line(sigvec,ncsvec);
      signoise_subplot(sigvec, ncsvec, m,b, ['S' num2str(i) ' Signal & Noise'])
   end

   subplot(2,2,4)
   ncsvec = get_upper(ncs_avg);
   [m, b] = fit_orthog_line(sigvec,ncsvec);
   signoise_subplot(sigvec, ncsvec, m,b, 'Avg. Signal & Noise')


end

function [] = signoise_subplot(sigvec,ncsvec,m,b,name)
   plot(sigvec, ncsvec, 'o'); hold on
   plot([-1,1], [-m+b, m+b], '--k')   
   grid on
   title(name)
   xlabel('Signal Corr.')
   ylabel('Noise Corr.')
   ylim([-0.5,1])
   legend({'Data', 'Fit line'}, 'Location', 'SouthEast')
end


function [] = plot_normed_sigcorr(sig_normed)

   figure()
   set(gcf,'Position', [100   517   458   427])
   scatter3(sig_normed(1,:), sig_normed(2,:), sig_normed(3,:));
   hold on
   xlim([0,1])
   ylim([0,1])
   zlim([0,1])

   [sphere_x,sphere_y,sphere_z] = sphere(100);
   color = ones(length(sphere_x));
   sphere_surf = surf(sphere_x,sphere_y,sphere_z, 'LineStyle', 'none', 'FaceAlpha', 0.3);

   title('Normed Signal Responses by Neuron')
   xlabel('S1 Response')
   ylabel('S2 Response')
   zlabel('S3 Response')
   
   set(gca, 'CameraPosition', [6.9846    2.3716    5.9265])
end

function [] = plot_total_corr_ranks(ranks1,ranks2,ranks3)

   figure()
   set(gcf,'Position', [95   274   446   670])

   % Permutation matrix which takes ranks to inds (i.e. st P*ranks' = inds) is P = inds'*ranks.
   P1 = ranks2' == ranks1;
   P2 = ranks3' == ranks1;

   [rows1, ~] = find(P1);
   [rows2, ~] = find(P2);

   orders = [(1:11)', rows1, rows2];

   plot(orders', '-o');
   grid on
   xlim([0.5,3.5])
   ylim([0,12])
   xticks([1,2,3])

   xlabel('Stimulus #')
   ylabel('Rank Order')
   title('Total Correlation by Neuron')

end