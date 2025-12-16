

% fig 4a: advantage of non-isotropic noise
    % define error function
    rng(2);
    err_fun = @(A, B) mean((A - B).^2, 'all') / mean(B.^2, 'all');

    % simulate with RRR
    % maximum variance 
    M_sigma = linspace(1,500,6);

    % run simulation
    n_sim = 50;
    ny = 50; nx=50;
    err = nan(n_sim,numel(M_sigma),2);
    for ii = 1:numel(M_sigma)
        % generate ground truth covariance
        Sigma = diag([ones(1,ny-1) M_sigma(ii)]); % variance all ones but the last 

        for isim = 1:n_sim
            [X,Y,U,V,ops] = simu_RRR(struct('Sigma',Sigma,'T',1000,'nx',nx,'ny',ny,'signse',50));

            % estimate, cross validation
            dW = Inf;
            w0 = svd_RRR(X, Y, ops.rnk);
            W_prev = w0;
            while dW > 1e-6
                cov_Y_res = cov(Y - X*W_prev);
                w0_noniso = svd_RRR_noniso(X, Y, ops.rnk,cov_Y_res);
                dW = norm(w0_noniso - W_prev,'fro');
                W_prev = w0_noniso;
                % fprintf('dW = %.6f\n',dW);
            end
            err(isim,ii,:) = cellfun(@(x) err_fun(x,U*V'),{w0,w0_noniso});
        end
    end

    % plot
    close all;
    Colors = colormap('lines');
    ax = subplot(1,1,1,'NextPlot','add');
    for ii = 1:2
        plot_multiple_lines(err(:,:,ii),ax,'x',M_sigma,'base_color',Colors(ii,:),'to_plt',{'err'});
    end
    xlabel('max(\sigma^2)','Interpreter','tex');
    ylabel('estimation error');
    t = text_legend(ax,{'standard RRR','full-covariance RRR'},[],1);
    set(findall(gcf,'-property','FontSize'),'FontSize',10);
    set(ax,'FontSize',9,'XLim',[min(M_sigma) max(M_sigma)]);
    set(gcf,'Position',[0 0 3 2.2]*72);

    saveas(gcf,'Fig4A','png');

    return



% fig 3a: advantage of regularization 

    % amount of data, estimation error w/ or w/o regularization
    n_sim = 50;  
    nx = 50;ny = 50;
    NSamples = round(50*2.^(1:.5:4));
    % Lambdas = 0:0.1:2;
    Lambdas = [0:10 20 100 1e3];
    err_fun = @(A, B) mean((A - B).^2, 'all') / mean(B.^2, 'all') ;

    U = [];V = [];
    err = nan(n_sim,numel(Lambdas),numel(NSamples));
    for isim = 1:n_sim
        [X,Y,U,V,ops] = simu_RRR(struct('rnk',2,'nx',nx,'ny',ny,'signse',50,'T',max(NSamples),'U',U,'V',V));
        w_true = U*V';
        for jj = 1:numel(Lambdas)
            w_hat = arrayfun(@(n_sample) svd_RRR(X(1:n_sample,:), Y(1:n_sample,:), ops.rnk,Lambdas(jj)),NSamples,'uni',0);
            err(isim,jj,:) = cellfun(@(w_) err_fun(w_,w_true), w_hat);
        end
    end


    % identify best lambda
    tmp = squeeze(mean(err,1));
    [~,I] = min(tmp(:,2));
    % plotting 
    close all;
    ax = subplot(1,1,1,'NextPlot','add');
    % plot best & no lambda
    plot_multiple_lines(permute(err(:,[1 I],:),[1 3 2]),ax(1),'to_plt',{'err'},'x',NSamples);
    text_legend(ax(1),{'standard RRR','Ridge-RRR'},[],2);ylabel(ax(1),'estimation error');xlabel(ax(1),'# samples');
    % plot all lambdas
    set(findall(gcf,'-property','FontSize'),'FontSize',10);
    set(ax,'FontSize',9,'XLim',[min(NSamples) max(NSamples)]);
    set(gcf,'Position',[0 0 3 2.2]*72);
    saveas(gcf,'Fig3A','png');
    return

% fig 6: output alignment

    T = 500;
    nx = 20;
    ny = 15;
    X = randn(T,nx);

    rnk = 3;
    U = randn(nx,rnk);
    V = randn(ny,rnk);
    Y_comm_only = X*U*V' + randn(T,ny) * 0;  % Y is aligned with X in the direction of U*V'
    

    % use true W
    W = U * V';

    % just do covariance
    basis = orth(randn(ny,ny));
    tmp = linspace(0, 1, ny);  % normalized time from 0 to 1
    start_val = 2;
    end_val = 1;
    s_y = end_val + (start_val - end_val) * exp(-3 * tmp);  % adjust -3 for steepness
    cov_Y = basis' * diag(s_y) * basis;

    
    
    % plot
    close all;
    r = 3; c = 3; ax = arrayfun(@(ii) subplot(r,c,ii,'NextPlot','add'),1:(r*c));
    set(groot, 'defaultLineMarkerSize', 3);
    for ii = 1:3
        switch ii
        case 1 % align with output
            Y = Y_comm_only + randn(T,ny) * 10; % add small noise
            cov_XW = basis' * diag([s_y(1:rnk) zeros(1,ny-rnk)]) * basis;
        case 2 % misalign with output
            tmp = randn(ny,ny-rnk);
            tmp = tmp - Y_comm_only' * (pinv(Y_comm_only') * tmp); % remove the component in the direction of U*V'
            tmp = tmp ./ vecnorm(tmp);
            Y = Y_comm_only + 1e3*randn(T,ny-rnk) * tmp';  % add noise in the direction of the 2nd principal component
            cov_XW = basis' * diag([zeros(1,ny-rnk) s_y(ny-rnk+1:ny)]) * basis;
        case 3 % random
            Y = Y_comm_only + randn(T,ny) * 1e3; % add large noise
            cov_XW = basis' * diag(s_y .* rand(1,ny)) * basis;
        end

        % plot output alignment
        [upop,spopvec] = svd(cov_Y,'vector'); 
        spopcum = cumsum(spopvec);
        spopvec_nrm = spopvec/sum(spopvec);
        spopcum_nrm = cumsum(spopvec_nrm);
        % Project communication covariance onto PCs
        cov_predicted = cov_XW; 
        scomvec = diag(upop'*cov_predicted*upop);
        scomvec_nrm = scomvec./sum(scomvec);
        scomcum_nrm = cumsum(scomvec)/sum(scomvec);
        nd = length(spopvec);
        muscom = mean(scomcum_nrm);
        muspop = mean(spopcum_nrm);
        clrs = get(gca,'colororder');
        commfrac = sum(scomvec)/sum(spopvec);

        % compute alignment index
            muspop = mean(spopcum_nrm);  % mean of PCA cumulative eigenvalues
            muscom = mean(scomcum_nrm);  % mean of communication cumulative eigenvalues

            alignment_raw = muscom-muspop;  % raw alignment score

            % compute communication fraction
            commfrac = sum(scomvec)/sum(spopvec);

            % compute max possible alignment
            totcom = sum(scomvec);
            ind = find(spopcum>totcom+1e-10,1);  % find how many dimensions we'd need for maximally aligned communication
            scommax = spopvec;
            scommax(ind:end) = 0;
            scommax(ind) = totcom-sum(scommax); % fix up final bin
            scommax_cum = cumsum(scommax)/sum(scommax); % compute cumulative
            a_max = mean(scommax_cum)-muspop; % max possible alignment score

            % compute min possible alignment (by flipping order of pop eigenvalues)
            spopvec_rev = flipud(spopvec);
            spopcum_rev = cumsum(spopvec_rev);
            ind = find(spopcum_rev>totcom+1e-10,1);  % find how many dimensions we'd need for minimally aligned communication
            scommin = spopvec_rev;
            scommin(ind:end) = 0;
            scommin(ind) = totcom-sum(scommin); % fix up final bin
            scommin = flipud(scommin); % flip back to normal ordering
            scommin_cum = cumsum(scommin)/sum(scommin); % compute cumulative
            a_min = mean(scommin_cum)-muspop; % min possible alignment score
            % rescale alignment scores above and below zero
            if alignment_raw>0
                a_out = alignment_raw/a_max;
            else
                a_out = alignment_raw/abs(a_min);
            end

        axes(ax(sub2ind([c r],1,ii)));
        plot(1:nd,spopvec_nrm,'o-',1:nd,scomvec/sum(spopvec),'o-');
        set(gca,'ylim',[0 spopvec_nrm(1)*1.1],'XLim',[0.5 ny+0.5]);

        axes(ax(sub2ind([c r],2,ii)));
        plot(1:nd,spopcum_nrm,'-o',1:nd,scomcum_nrm,'o-'); box off;
        set(gca,'ylim',[0 1],'YTick',[0 1],'XLim',[.5 ny+0.5]);

        axes(ax(sub2ind([c r],3,ii)));
        plot(commfrac,a_out,'k*'); 
        plot([0 1],[0 0],'k--', [0.5 0.5], [-1 1], 'k--','linewidth',1); 
        set(gca,'ylim',[-1.01 1.01], 'xtick',-1:.5:1,'xlim',[0 1],'YTick',[-1 0 1]);
        axis square;
        
        box off;
    end
    % % find line plot in ax(1)
    % h = findall(ax(1),'Type','Line');
    
    % axis label
    ylabel(ax(1),'normalized variance'); 
    ylabel(ax(2),{'normalized','cumulative variance'});
    ylabel(ax(3),{'output alignment','index'});
    xlabel(ax(9),'communication fraction');
    xlabel(ax(7),'PC dimension');
    xlabel(ax(8),'PC dimension');
    set(ax(1:6),'XTick',[]);
    % get all labels
    set(findall(gcf,'-property','FontSize'),'FontSize',9);
    set(ax,'FontSize',8);
    % find all axis labels
    set(findall(gcf,'-property','FontName'),'FontName','Arial');


    set(gcf,'Position',[0 0 6 5] * 72);
    saveas(gcf,'fig6','png');

    return


% figure 5: input alignment

    T = 1000; % number of time points
    nx = 20;
    ny = 15;
    r = 3; % rank of the communication subspace

   
    % create exponential decay eigenvalues
    basis = orth(randn(nx,nx)); % random orthogonal basis
    n_samples = 20;
    tmp = linspace(0, 1, n_samples);  % normalized time from 0 to 1
    start_val = 2;
    end_val = 1;
    s = end_val + (start_val - end_val) * exp(-3 * tmp);  % adjust -3 for steepness
    % generate input with fixed eigenvalues
    cov_X = basis * diag(s) * basis'; % covariance of input
    X = mvnrnd(zeros(nx,1),cov_X,T); % generate data

    % input main axes
    [Upca,Spcavec] = svd(cov_X,"vector"); % input axes


    Ws = {Upca(:,1:r), ...
        Upca(:,end-r+1:end), ...
        randn(nx,r) }; % random weights
    
    close all;
    ax = arrayfun(@(ii) subplot(1,3,ii,'NextPlot','add') ,1:3);
    for ii = 1:3
        % Compute SVD of weights
        W = Ws{ii}; % pick weight vector
        W = W / norm(W); % normalize weight (magnitude does not affect alignment)
        [Uw,Sw,Vw]= svd(W);
        Swvec = diag(Sw);
        var_comm = diag(Upca'*W*W'*Upca); % sum(Sw' * ((Uw'*Upca).^2) .* Spcavec',1)


        a_in = alignment_input(X, W , r, cov_X);

        % --------  Make plot --------------
        axes(ax(ii));
        yyaxis(ax(ii),'left');plot(Spcavec,'o-','MarkerSize',3);  %ylabel('input variance');
        yyaxis(ax(ii),'right');plot(var_comm,'o-','MarkerSize',3); %ylabel('communicated variance');
        ax(ii).YAxis(1).Limits = [-.01 max(Spcavec)*1.15];
        ax(ii).YAxis(2).Limits = [-.01 max(var_comm)*1.15];
        set(ax(ii),'XLim',[0.5 nx+0.5]);
        set(ax(ii).YAxis(1),'Limits',[0 max(Spcavec)*1.15],'TickValues',[0 1 2]);
        set(ax(ii).YAxis(2),'Limits',[-0.01 1.2],'TickValues',[0 1]);
        t = text(ax(ii).XLim(2),ax(ii).YAxis(2).Limits(2),...
            sprintf('input alignment = %.2f', a_in),...
            'HorizontalAlignment','right','VerticalAlignment','top',...
            'FontSize',8,'FontWeight','normal');
        
    end
    set(ax,'FontSize',8);
    set(ax(1).YAxis(1).Label,'String','input variance','FontSize',9);
    set(ax(3).YAxis(2).Label,'String','variance communicated','FontSize',9);
    xlabel(ax(2),'input PC dimension','FontSize',9);
    set(gcf,'Position',[0 0 7 1.8]*72);

    saveas(gcf,'fig5','png');

    return


