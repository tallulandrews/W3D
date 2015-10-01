# Modularize this stuff more sensibly
#  Plotting Functions
bg__dropout_plot_base <- function (norm, weights = 1, xlim = NA) {
	require("RColorBrewer")
	
	gene_info = calc_variables(norm,weights);

        xes = log(gene_info$s)/log(10);
        put_in_order = order(xes);
        fancy <- densCols(xes, gene_info$p, colramp=colorRampPalette(c("black","white")))
        dens <- col2rgb(fancy)[1,]+1L
        colours <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F",
                                    "#FCFF00", "#FF9400", "#FF3100"))(256)
        dens.col = colours[dens]

        par(fg="black")
	if (!is.na(xlim)) {
        	plot(xes,gene_info$p, main=name, ylab="Dropout Proportion", xlab="log(expression)", col = dens.col,pch=16, xlim=xlim, ylim=c(0,1))
	} else {
        	plot(xes,gene_info$p, main=name, ylab="Dropout Proportion", xlab="log(expression)", col = dens.col,pch=16)
	}
	return(list(P=gene_info$p, S=gene_info$s, xes=xes, data=norm, weights=weights, order=put_in_order));
}

bg__add_model_to_plot <- function(fitted_model, base_plot, lty=1, lwd=1, col="black",legend_loc = "topright") {
	lines(base_plot$xes[base_plot$order],fitted_model$predictions[base_plot$order],lty=lty,lwd=lwd,col=col);
	par(fg=col)
        this_loc = legend(legend_loc, fitted_model$model, box.lty=lty, box.lwd=lwd, xjust=1)
	par(fg="black")
	return(this_loc)
}

bg__highlight_genes <- function (base_plot, genes, colour="purple", pch=16) {
	if(!is.numeric(genes)) {
		genes = match(genes, rownames(base_plot));
		nomatch = sum(is.na(genes));
		if (nomatch > 0) {warning(paste(nomatch, " genes could not be matched to data, they will not be highlighted."));}
		genes = genes[!is.na(genes)];
	}
	points(base_plot$xes[genes],base_plot$p[genes],col=colour, pch=pch)
}

bg__expression_heatmap <- function (genes, data, cell_labels=NA, gene_labels=NA) { # Should weighting be taken into account for heatmap clustering?
	require("RColorBrewer")
	require("gplots")
	if(!is.numeric(genes)) {
		genes = match(genes, rownames(base_plot));
		nomatch = sum(is.na(genes));
		if (nomatch > 0) {warning(paste(nomatch, " genes could not be matched to data, they will not be included in the heatmap."));}
		genes = genes[!is.na(genes)];
	}
	# Plot heatmap of expression
	heatcolours <- brewer.pal(11,"RdBu")
	heat_data = as.matrix(data[genes,])
	heat_data = log(heat_data+1)/log(2);
	ColColors = rep("white", times=length(heat_data[1,]))
	RowColors = rep("white", times=length(heat_data[1,]))
	if (!is.na(cell_labels)) {
		colours = as.factor(cell_labels)
		palette = brewer.pal(max(3,length(unique(cell_labels))), "Set3");
		ColColors = palette[colours];	
	} 
	if (!is.na(gene_labels)) {
		# lowest factor level = grey (so 0-1 is striking)
		colours = as.factor(gene_labels)
		palette = c("grey75",brewer.pal(max(3,length(unique(gene_labels))), "Set1"));
		RowColors = palette[colours];
	}
	heatmap.2(heat_data, ColSideColors = ColColors, RowSideColors = RowColors, col=heatcolours, scale="row",symbreaks=T, trace="none", dendrogram="column")
}

# Model-fitting/manipulation Functions
bg__calc_variables <- function(norm, weights = 1) {
        # Calc variables
	if (prod(dim(weights) == dim(norm))) {
		p = rowZero_wgt(norm,weights)/rowSums(weights);
        	s = rowMeans_wgt(norm,weights);
		s_stderr = sqrt(rowVar_wgt(norm,weights))/sqrt(rowSums(weights));
		p_stderr = sqrt(p_obs*(1-p_obs)/rowSums(weights));
	} else {
		print("Weights not provided or not same dimension as expression matrix. Using unweighted version.")
	        p = apply(norm,1,function(x){y = x[!is.na(x)]; sum(y==0)/length(y)});
	        s = rowMeans(norm, na.rm=T);
		s_stderr = unlist(apply(norm,1,sd));
		s_stderr = unlist(apply(norm,1,sd))/sqrt(length(norm[1,]));
		p_stderr = sqrt(p_obs*(1-p_obs)/length(norm[1,]));
	}
	names(s) = rownames(norm);
	names(p) = rownames(norm);
	return(list(s = s, p = p, s_stderr = s_stderr, p_stderr = p_stderr))
}

bg__invert_MM <- function (K, p) {K*(1-p)/(p)}
bg__horizontal_residuals_MM_log10 <- function (K, p, s) {log(s)/log(10) - log(bg__invert_MM(K,p))/log(10)}
bg__num.zero <- function(x){sum(x==0)}
bg__fit_MM <- function (p,s) {
        fit = nls(p ~ 1-(s/((krt+s))),data.frame(s=s),start=list(krt=3))
	K_glm = glm(p ~ offset(-1*log(s)), family="binomial")
	Kerr = summary(K_glm)$coeff[1,2];
	Kcoeff = summary(K_glm)$coeff[1,1];
	Kerr = exp(Kcoeff+Kerr)-exp(Kcoeff)
        predicted = fitted(fit)
        krt=summary(fit)$parameters[1,1]
	return(list(K=krt,Kerr=Kerr,predictions=predicted, model=c("MMenton",paste("Krt =",round(krt,digits=3))),SSr=round(sum((residuals(doubleXfit))^2)),SAr=round(sum(abs(residuals(doubleXfit))))))
}
bg__fit_logistic <- function(p,s) {
        logistic = glm(p~log(s),family="binomial")
        predlog = fitted(logistic)
	return(list(predictions=predlog, model=c( "Logistic", paste("Intercept =",round(logistic$coeff[1],digits=3)),paste("Coeff =",round(logistic$coeff[2],digits=3))),SSr=round(sum((fitted(logistic)-p)^2))),SAr=round(sum(abs(fitted(logistic)-p))));
}
bg__fit_ZIFA <- function(p,s) {
	doubleXfit = nls(p ~ exp(-lambda*s*s),data.frame(s=s),start=list(lambda=0), algorithm="port", lower=list(lambda=0));
	preddoubleX = fitted(doubleXfit);
	lambda=summary(doubleXfit)$parameters[1,1];
	return(list(predictions=preddoubleX, model=c("p ~ e^(-lambda*S^2)",paste("lambda =",signif(lambda,digits=2))),SSr = round(sum((residuals(doubleXfit))^2)),SAr = round(sum(abs(residuals(doubleXfit))))));
}

# Normalization Functions
bg__UQ <- function(x){quantile(x[x>0],0.75)};
bg__filter_genes <- function(data) {
        # get rid of genes with 0 expression
        filter <- apply(data, 1, function(x) length(x[x>5])>=2);
        data = data[filter,];
	return(data);
}

bg__filter_cells <- function(data,labels=NA, suppress.plot=FALSE) {
	cell_zero = apply(data,2, bg__num.zero)/length(data[,1]);
	mu = mean(cell_zero);
	sigma = sd(cell_zero);
	# Deal with bi-modal
	if (sum(cell_zero > mu-sigma & cell_zero < mu+sigma) < 0.5) { # should be 0.68 theoretically
		mu = mean(cell_zero[cell_zero < median(cell_zero)]);
		sigma = sd(cell_zero[cell_zero < median(cell_zero)]);
	}
	low_quality = p.adjust(pnorm((cell_zero-mu)/sigma, lower.tail=F), method="fdr") < 0.05;
	if (sum(low_quality) > 0) {
		data = data[,!low_quality];
		cell_zero = cell_zero[!low_quality];
		if (!is.na(labels)) {labels = labels[!low_quality]}
	}
	if (!suppress.plot) {
		hist(cell_zero, col="grey75", xlab="Number of zeros (per cell)", main="", prob=TRUE)
		curve(dnorm(x,mean=mu, sd=sigma), add=TRUE)
		if (sum(low_quality) > 0) {
			abline(v=min(cell_zero[low_quality]), col="red")
		}
	}
	return(list(data = data, labels = labels));
}

bg__normalize <- function(data) {
	# Combine UQ and detection rate adjusted normalization 
	# Stephanie Hick, Mingziang Teng, Rafael A Irizarry "On the widespread and critical impact of systematic single-cell RNA-Seq data" http://dx.doi.org/10.1101/025528 
	uq = unlist(apply(data,2,bg__UQ));
	normfactor = (uq/median(uq)) * (median(cell_zero)/cell_zero); 
	data = t(t(data)/normfactor);
	return(data);
}

# DE Genes functions

# Use the fact that errors of proportions are well define by converting S to proportion detected equivalents?
obsolete__test_DE_P_equiv <- function (norm, weights=1, fit=NA) {
	gene_info = bg__calc_variables(norm, weights);
	if (is.na(fit)) {
		fit = bg__fit_MM(gene_info$p, gene_info$s);
	}
	p_obs = gene_info$p;
	N = length(norm[1,]);
	p_err = gene_info$p_stderr;
	S_mean = gene_info$s
	S_err = gene_info$s_stderr
	K_err = fit$Kerr;
	p_equiv = fit$K/(fit$K+S_mean);
	propagated_err_p_equiv = p_equiv*sqrt(((S_err+K_err)/(S_mean+fit$K))^2+(K_err/fit$K)^2)
	Z = (p_equiv - p_obs)/sqrt(p_err^2+propagated_err_p_equiv^2); # low = shifted right, high = shifted left
	pval = pnorm(Z, lower.tail=T)
	effect_size = p_obs/p_equiv;
	return(list(pval = pval, fold_change = effect_size))
}

# Use the fact that S as a function of P is more stable to noise for the main part of the curve
bg__test_DE_S_equiv <- function (norm, weights=1, fit=NA, method="propagate") {
	gene_info = bg__calc_variables(norm, weights);
	if (is.na(fit)) {
		fit = bg__fit_MM(gene_info$p, gene_info$s);
	}
	p_obs = gene_info$p;
	N = length(norm[1,]);
	p_err = gene_info$p_stderr;
	S_mean = gene_info$s
	S_err = gene_info$s_stderr
	K_err = fit$Kerr;
	S_equiv = bg__invert_model(fit$K,p_obs);

	## Monte Carlo method to estimate error around S_equiv ##
	MC_err <- function (p_base) {
		p_rand = rnorm(10000, mean = p_base, sd = p_err);
		p_rand = p_rand[p_rand > 0 & p_rand < 1]
		K_rand = rnorm(length(p_rand),fit$K,sd = K_err);
		K_rand[K_rand < 1] = 1;
		S_equiv_rand = bg__invert_model(K_rand, p_rand)
		sd(S_equiv_rand)
	}
	if (method == "MC") {
		S_equiv_err = unlist(lapply(p_obs,MC_err))
	} else {
		S_equiv_err = S_equiv*sqrt(2*(p_err/p_obs)^2+(K_err/fit$K)^2);
	}

	Z = (S_equiv - S_mean)/sqrt(S_err^2+S_equiv_err^2); # low = shifted right, high = shifted left
	pval = pnorm(Z, lower.tail=T)*2
	effect_size = (S_mean-S_equiv)/S_equiv;
	return(list(pval = pval, effect = effect_size))
}

bg__get_extreme_residuals <- function (norm,weights, v_threshold=c(0.05,0.95), fdr_threshold = 0.1, direction="right", suppress.plot = FALSE) {
	gene_info = bg__calc_variables(norm, weights);
	if (is.na(fit)) {
		fit = bg__fit_MM(gene_info$p, gene_info$s);
	}
	res = bg__horizontal_residuals_MM_log10(fit$K, gene_info$p, gene_info$s)
	res = res[gene_info$p < max(v_threshold) & gene_info$p > min(v_threshold)]
	mu = mean(res); sigma = sd(res);
	# deal with potential bi-modality
	if (sum(res > mu-sigma & res < mu+sigma) < 0.5) { # should be 0.68 theoretically
		mu = mean(cell_zero[cell_zero > quantile(cell_zero,0.33)]);
		sigma = sd(cell_zero[cell_zero > quantile(cell_zero,0.33)]);
	}

	if (direction="right") {
		pval =pnorm((res-mu)/sigma, lower.tail=F)
	} else {
		pval = pnorm((res-mu)/sigma, lower.tail=T)
	}
	qval = p.adjust(pval, method="fdr");
	sig = qval < fdr_threshold;

	# Plot fitted normal curve
	if (!suppress.plot) {
		hist(res, col="grey75", xlab="horizontal residuals", main="", prob=TRUE)
		curve(dnorm(x,mean=mu, sd=sigma), add=TRUE);
		if (direction="right" & sum(sig) > 0) {
			abline(v=min(res[sig]), col="red");
		} else {
			abline(v=max(res[sig]), col="red");
		}
	}
	return(names(pval)[sig]);
}
##### Assembled Analysis Chunks ####
W3D_Clean_Data <- function(data, labels = NA, suppress.plot=FALSE) {
	data = bg__filter_genes(data);
	data = bg__filter_cells(data, labels, suppress.plot = suppress.plot);
	norm = bg__normalize(data$data);
	return(list(data = norm, labels=labels));
}

W3D_Dropout_Models <- function(data_list, weights = 1, xlim=c(-1.5,6)) {
	BasePlot = bg__dropout_plot_base(data_list$data, weights = weights, xlim = xlim);
	MM = bg__fit_MM(BasePlot$P, BasePlot$S);
	SCDE = bg__fit_logistic(BasePlot$P, BasePlot$S);
	ZIFA = bg__fit_ZIFA(BasePlot$P, BasePlot$S);
	sizeloc = bg__add_model_to_plot(MM, BasePlot, lty=1, lwd=2.5, col="black",legend_loc = "topright");
	sizeloc = bg__add_model_to_plot(SCDE, BasePlot, lty=2, lwd=2.5, col="grey35",legend_loc = c(sizeloc$rect$left+sizeloc$rect$w,sizeloc$rect$top-sizeloc$rect$h-0.05));
	sizeloc = bg__add_model_to_plot(ZIFA, BasePlot, lty=3, lwd=2.5, col="grey75",legend_loc = c(sizeloc$rect$left+sizeloc$rect$w,sizeloc$rect$top-sizeloc$rect$h-0.05));
	return(list(MMfit = MM, LogiFit = SCDE, ExpoFit = ZIFA));
}

W3D_Differential_Expression <- function(data_list, weights, knownDEgenes=NA, xlim=c(-1.5,6), method="propagate") {
	BasePlot = bg__dropout_plot_base(data_list$data, weights = weights, xlim = xlim);
	MM = bg__fit_MM(BasePlot$P, BasePlot$S);
	sizeloc = bg__add_model_to_plot(MM, BasePlot, lty=1, lwd=2.5, col="black",legend_loc = "topright");
	DEgenes = bg__test_DE_S_equiv(data_list$data, weights=weights, fit=MM, method=method);
	bg__highlight_genes(DEgenes, BasePlot);
	bg__expression_heatmap(genes, data_list$data, cell_labels=data_list$labels, gene_labels=knownDEgenes);
	return(DEgenes)
}

W3D_Get_Extremes <- function(data_list, weights) {
	BasePlot = bg__dropout_plot_base(data_list$data, weights = weights, xlim = c(-1.5,6));
	MM = bg__fit_MM(BasePlot$P, BasePlot$S);
	shifted_right = bg__get_extreme_residuals(data_list$data,weights, v_threshold=c(0.05,0.95), fdr_threshold = 0.1, direction="right", suppress.plot=TRUE)
	shifted_left  = bg__get_extreme_residuals(data_list$data,weights, v_threshold=c(0.05,0.95), fdr_threshold = 0.1, direction="left",  suppress.plot=TRUE)
	bg__highlight_genes(shifted_right, BasePlot, col="orange");
	bg__highlight_genes(shifted_left,  BasePlot, col="purple");
	return(list(left=shifted_left,right=shifted_right));
}
