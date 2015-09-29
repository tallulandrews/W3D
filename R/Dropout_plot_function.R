
source("Weight_Dropouts.R")
require("RColorBrewer")
## Short functions ##
invert_model <- function (K, p) {K*(1-p)/(p)}
horizontal_residuals_log10 <- function (K, p, s) {log(s)/log(10) - log(invert_model(K,p))/log(10)}
num.zero <- function(x){sum(x==0)}
UQ <- function(x){quantile(x[x>0],0.75)};

## Added 9-10 Sept 2015 ##

get_extreme_residuals <- function (dropout_obj, threshold = 0.1, labels=NA) {
	res = horizontal_residuals_log10(dropout_obj$K, dropout_obj$P, dropout_obj$S)
	res = res[dropout_obj$P < 0.95 & dropout_obj$P > 0.05]
	hist(res)
	mu = mean(res); sigma = sd(res);
	pval =pnorm((res-mu)/sigma, lower.tail=F)
	qval = padjust(pval, method="fdr");
	sig = names(pval)[qval < threshold];
	sig.rows = match(sig,rownames(dropout_obj$data));
	# Plot DE genes on dropout plot
	dropout_plot(data,"DE Genes", c(-1.5,6));
	points(log(dropout_obj$S[sig.rows])/log(10), dropout_obj$P[sig.rows], col="black", pch=16)
	
	table = data.frame( gene = sig)
	for (tag in unique(labels)) {
		table = cbind(table,rowMeans(dropout_obj$data[sig.rows,labels==tag]));
	}
	table = cbind(table, qval);
	colnames(table) <- c("gene", unique(labels),qval);

	# Plot heatmap of expression
	heatcolours <- brewer.pal(11,"RdBu")
	if (is.na(labels)) {
		X11()
		heatmap(log(as.matrix(data[sig.rows,])+1), col=heatcolours)
	} else {
		colours = as.factor(labels)
		palette = topo.colours(n=length(unique((colours))));
		X11()
		heatmap(log(as.matrix(data[sig.rows,])+1), ColSideCol = palette[colours], col=heatcolours)
	}
	return(table)

}
get_DE_genes <- function(dropout_obj, mt_method="bon",threshold = 0.05, labels = NA) {
	require("gplots")
	DEgenes = test_DE_S_equiv(dropout_obj);
	pval = DEgenes$pval;
	pval[is.na(pval)] = 1;
	qval = p.adjust(pval, method=mt_method);
	sig = names(pval)[qval < threshold];
#	sig = names(pval)[pval < threshold/length(pval)];
	sig.rows = match(sig,rownames(dropout_obj$data));
	# Plot DE genes on dropout plot
#	png("Monocle_Dropout_Plot_LabMeeting2.png", width=6,height=6, units="in", res=300)
	blank = dropout_plot(dropout_obj$data,"DE Genes", c(-1.5,6));
	points(log(dropout_obj$S[sig.rows])/log(10), dropout_obj$P[sig.rows], col="purple", pch=16)
#	dev.off()

	table = data.frame( gene = sig)
	for (tag in unique(labels)) {
		table = cbind(table,rowMeans(dropout_obj$data[sig.rows,labels==tag]));
	}
	table = cbind(table, qval[sig.rows]);
	colnames(table) <- c("gene", unique(labels),"qval");
#	X11()
#	matplot(t(table[,2:5]), type="l", log="y", ylab="median expression", xlab="time (days)")
	# Plot heatmap of expression
	heatcolours <- brewer.pal(11,"RdBu")
	heat_data = as.matrix(dropout_obj$data[sig.rows,])
	heat_data = log(heat_data+1)
	if (is.na(labels)) {
		X11()
		heatmap.2(heat_data, col=heatcolours, scale="row",symbreaks=T, trace="none")
	} else {
		colours = as.factor(labels)
		palette = brewer.pal(max(3,length(unique(labels))), "Accent");
		X11()
#		png("Monocle_Dropout_Plot_LabMeeting3.png", width=6,height=6, units="in", res=300)
		heatmap.2(heat_data, ColSideCol = palette[colours], col=heatcolours, scale="row",symbreaks=T, trace="none")
#		dev.off()
	}
	return(table)
}
get_DE_genes_wgt <- function(dropout_obj, mt_method="bon", threshold = 0.05, labels = NA) {
	require("gplots")
	DEgenes = test_DE_S_equiv(dropout_obj);
	pval = DEgenes$pval;
	pval[is.na(pval)] = 1;
	qval = p.adjust(pval, method=mt_method);
	sig = names(pval)[qval < threshold];
#	sig = names(pval)[pval < threshold/length(pval)];
	sig.rows = match(sig,rownames(dropout_obj$data));
	# Plot DE genes on dropout plot
#	png("Monocle_Dropout_Plot_LabMeeting2.png", width=6,height=6, units="in", res=300)
	blank = dropout_plot_wgt(dropout_obj$data,"DE Genes", c(-1.5,6));
	points(log(dropout_obj$S[sig.rows])/log(10), dropout_obj$P[sig.rows], col="purple", pch=16)
#	dev.off()

	table = data.frame( gene = sig)
	for (tag in unique(labels)) {
		table = cbind(table,rowMeans(dropout_obj$data[sig.rows,labels==tag]));
	}
	table = cbind(table, qval[sig.rows]);
	colnames(table) <- c("gene", unique(labels),"qval");
#	X11()
#	matplot(t(table[,2:5]), type="l", log="y", ylab="median expression", xlab="time (days)")
	# Plot heatmap of expression
	heatcolours <- brewer.pal(11,"RdBu")
	heat_data = as.matrix(dropout_obj$data[sig.rows,])
	heat_data = log(heat_data+1)
	if (is.na(labels)) {
		X11()
		heatmap.2(heat_data, col=heatcolours, scale="row",symbreaks=T,trace="none")
	} else {
		colours = as.factor(labels)
		palette = brewer.pal(max(3,length(unique(labels))), "Accent");
		X11()
#		png("Monocle_Dropout_Plot_LabMeeting3.png", width=6,height=6, units="in", res=300)
		heatmap.2(heat_data, ColSideCol = palette[colours], col=heatcolours, scale="row",symbreaks=T, trace="none")
#		dev.off()
	}
	return(table)
}



# Use the fact that errors of proportions are well define by converting S to proportion detected equivalents?
test_DE_P_equiv <- function (dropout_obj) {
	p_obs = dropout_obj$P;
	N = length(dropout_obj$data[1,]);
	sd_p = sqrt(p_obs*(1-p_obs));
	S_mean = rowMeans(dropout_obj$data)
	S_sd = unlist(apply(dropout_obj$data,1,sd))
	K_sd = dropout_obj$Kerr;
	p_equiv = dropout_obj$K/(dropout_obj$K+S_mean);
	propagated_err_p_equiv = p_equiv*sqrt(((S_sd+K_sd)/(S_mean+dropout_obj$K))^2+(K_sd/dropout_obj$K)^2)
	Z = (p_equiv - p_obs)/sqrt(sd_p^2/N+propagated_err_p_equiv^2/N); # low = shifted right, high = shifted left
#	S_pred = invert_model(dropout_obj$K,dropout_obj$P)
#	Z = (S_mean-S_pred)/(S_sd/sqrt(length(dropout_obj$data[1,])))
	pval = pnorm(Z, lower.tail=T)
	effect_size = (p_obs-p_equiv)/p_equiv;
	return(list(pval = pval, effect = effect_size))
}

# Use the fact horizontal distances are less senstive to error by converting p to expression equivalents
# get NA if p = 0;
test_DE_S_equiv <- function (dropout_obj) {
	p_obs = dropout_obj$P;
	N = length(dropout_obj$data[1,]);
	sd_p = sqrt(p_obs*(1-p_obs));
	S_mean = rowMeans(dropout_obj$data)
	S_sd = unlist(apply(dropout_obj$data,1,sd));
	K_sd = dropout_obj$Kerr;
	S_equiv = invert_model(dropout_obj$K,dropout_obj$P);

	## Monte Carlo method to estimate error around S_equiv ##
#	MC_err <- function (p_base) {
#		p_rand = rnorm(10000, mean = p_base, sd = sqrt(p_base*(1-p_base)/N));
#		p_rand = p_rand[p_rand > 0 & p_rand < 1]
#		K_rand = rnorm(length(p_rand),dropout_obj$K,sd = K_sd);
#		K_rand[K_rand < 1] = 1;
#		S_equiv_rand = invert_model(K_rand, p_rand)
#		sd(S_equiv_rand)
#	}
#	S_equiv_err = unlist(lapply(p_obs,MC_err))

	propagated_err_S_equiv = S_equiv*sqrt(2*(sd_p/p_obs)^2+(K_sd/dropout_obj$K)^2);

#	biggest_err = unlist(apply(cbind(S_equiv_err,propagated_err_S_equiv),1,max));

	Z = (S_equiv - S_mean)/sqrt(S_sd^2/N+propagated_err_S_equiv^2/N); # low = shifted right, high = shifted left
	pval = pnorm(Z, lower.tail=T)*2
	effect_size = (S_mean-S_equiv)/S_equiv;
	return(list(pval = pval, effect = effect_size))
}

test_DE_simple <- function (dropout_obj) {
	p_obs = dropout_obj$P;
        N = length(dropout_obj$data[1,]);
        sd_p = sqrt(p_obs*(1-p_obs));
        S_mean = rowMeans(dropout_obj$data)
        S_sd = unlist(apply(dropout_obj$data,1,sd))
        K_sd = 0;
        S_equiv = invert_model(dropout_obj$K,dropout_obj$P);
	Z = (S_mean-S_equiv)/(S_sd/sqrt(N));
        pval = pnorm(Z, lower.tail=F)
        return(pval)
}

## -------------------  ##
filter_normalize <- function(data, labels = NA) {
        # get rid of genes with 0 expression
        filter <- apply(data, 1, function(x) length(x[x>5])>=2);
        data = data[filter,];
	cell_zero = apply(data,2, num.zero)/length(data[,1]);
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

	# Combine UQ and detection rate adjusted normalization 
	# Stephanie Hick, Mingziang Teng, Rafael A Irizarry "On the widespread and critical impact of systematic single-cell RNA-Seq data" http://dx.doi.org/10.1101/025528 
	uq = unlist(apply(data,2,UQ));
	normfactor = (uq/median(uq)) * (median(cell_zero)/cell_zero); 
	data = t(t(data)/normfactor);
	weights = Calc_dropout_weights(data);
	return(list(data = data, factors = normfactor, labels = labels, weights = weights));
}	

dropout_plot <- function(data,name,xlim, labels = NA) {
	require("RColorBrewer")
        # get rid of genes with 0 expression
# Change filtering 9 Sept 2015
#        filter <- apply(data, 1, function(x) length(x[x>5])>=2)
#        data = data[filter,]
	if (length(labels) > 1 | !is.na(labels)) {
		norm = filter_normalize(data, labels);
	} else {
		norm = filter_normalize(data);
	}


        # Calc variables
        p = apply(norm$data,1,function(x){y = x[!is.na(x)]; sum(y==0)/length(y)});
#       p = p/length(data[1,]);
        s = rowMeans(norm$data, na.rm=T);

        # Fit Michalis Menten
        fit = nls(p ~ 1-(s/((krt+s))),data.frame(s=s),start=list(krt=3))#, algorithm="port", lower=list(krt=0));
	K_glm = glm(p ~ offset(-1*log(s)), family="binomial")
	Kerr = summary(K_glm)$coeff[1,2];
#	Kcoeff = summary(K_glm)$coeff[1,1];
#	Kerr = exp(Kcoeff+Kerr)-exp(Kcoeff)
	Kerr = exp(Kerr)

        # Fit Hill
#       fit = nls(p ~ 1-(1/(1+(kpcr/s)^n)),data.frame(s=s),start=list(kpcr=5,n=1), algorithm="port", lower=list(kpcr=0,n=0));
        predicted = fitted(fit)
        summary(fit)$parameters
        krt=summary(fit)$parameters[1,1]

        #SCDE
        logistic = glm(p~log(s),family="binomial")
        predlog = fitted(logistic)

        xes = log(s)/log(10)
        put_in_order = order(xes);
        fancy <- densCols(xes, p, colramp=colorRampPalette(c("black","white")))
        dens <- col2rgb(fancy)[1,]+1L
        colours <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F",
                                    "#FCFF00", "#FF9400", "#FF3100"))(256)
        dens.col = colours[dens]

        par(fg="black")
        #plot(xes,p, main=name, ylab="Dropout Proportion", xlab="log(expression)", col = dens.col,pch=16)
        plot(xes,p, main=name, ylab="Dropout Proportion", xlab="log(expression)", col = dens.col,pch=16, xlim=xlim, ylim=c(0,1))
        # Fits
        lines(xes[put_in_order],predicted[put_in_order], col="black", lwd=3, lty=1)
        lines(xes[put_in_order],predlog[put_in_order], col="grey50", lwd=3,lty=3)
        par(fg="grey50")
        sizeloc = legend("topright", c( "Logistic", paste("Intercept =",round(logistic$coeff[1],digits=3)),paste("Coeff =",round(logistic$coeff[2],digits=3))),box.lty=3, box.lwd=3)
        par(fg="black")
        legend(sizeloc$rect$left+sizeloc$rect$w,sizeloc$rect$top-sizeloc$rect$h-0.05, c("MMenton",paste("Krt =",round(krt,digits=3))), box.lty=1, box.lwd=3, xjust=1)
        if (length(grep("^ERCC",rownames(norm$data))) > 0) {
                erccs = grep("^ERCC",rownames(norm$data));
                points(xes[erccs],p[erccs],col="black",cex=2)
        }
	return(list(K = krt,P = p, S = s, data = norm$data, labels = norm$labels, Kerr = Kerr));
}


dropout_plot_wgt <- function(data,name,xlim, labels = NA) {
	require("RColorBrewer")
        # get rid of genes with 0 expression
# Change filtering 9 Sept 2015
#        filter <- apply(data, 1, function(x) length(x[x>5])>=2)
#        data = data[filter,]
	if (length(labels) > 1 | !is.na(labels)) {
		norm = filter_normalize(data, labels);
	} else {
		norm = filter_normalize(data);
	}


        # Calc variables
#        p = apply(norm$data,1,function(x){y = x[!is.na(x)]; sum(y==0)/length(y)});
	p = rowZero_wgt(norm$data,norm$weights)/rowSums(norm$weights)
#        s = rowMeans(norm$data, na.rm=T);
        s = rowMeans_wgt(norm$data,norm$weights);

        # Fit Michalis Menten
        fit = nls(p ~ 1-(s/((krt+s))),data.frame(s=s),start=list(krt=3))#, algorithm="port", lower=list(krt=0));
#        fit = nls(p ~ 1-(s/((krt+s))),data.frame(s=s),start=list(krt=3), algorithm="port", lower=list(krt=0));
	K_glm = glm(p ~ offset(-1*log(s)), family="binomial")
	Kerr = summary(K_glm)$coeff[1,2];
#	Kcoeff = summary(K_glm)$coeff[1,1];
#	Kerr = exp(Kcoeff+Kerr)-exp(Kcoeff)
	Kerr = exp(Kerr)

        # Fit Hill
#       fit = nls(p ~ 1-(1/(1+(kpcr/s)^n)),data.frame(s=s),start=list(kpcr=5,n=1), algorithm="port", lower=list(kpcr=0,n=0));
        predicted = fitted(fit)
        summary(fit)$parameters
        krt=summary(fit)$parameters[1,1]

        #SCDE
        logistic = glm(p~log(s),family="binomial")
        predlog = fitted(logistic)

        xes = log(s)/log(10)
        put_in_order = order(xes);
        fancy <- densCols(xes, p, colramp=colorRampPalette(c("black","white")))
        dens <- col2rgb(fancy)[1,]+1L
        colours <-  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F",
                                    "#FCFF00", "#FF9400", "#FF3100"))(256)
        dens.col = colours[dens]

        par(fg="black")
        #plot(xes,p, main=name, ylab="Dropout Proportion", xlab="log(expression)", col = dens.col,pch=16)
        plot(xes,p, main=name, ylab="Dropout Proportion", xlab="log(expression)", col = dens.col,pch=16, xlim=xlim, ylim=c(0,1))
        # Fits
        lines(xes[put_in_order],predicted[put_in_order], col="black", lwd=3, lty=1)
        lines(xes[put_in_order],predlog[put_in_order], col="grey50", lwd=3,lty=3)
        par(fg="grey50")
        sizeloc = legend("topright", c( "Logistic", paste("Intercept =",round(logistic$coeff[1],digits=3)),paste("Coeff =",round(logistic$coeff[2],digits=3))),box.lty=3, box.lwd=3)
        par(fg="black")
        legend(sizeloc$rect$left+sizeloc$rect$w,sizeloc$rect$top-sizeloc$rect$h-0.05, c("MMenton",paste("Krt =",round(krt,digits=3))), box.lty=1, box.lwd=3, xjust=1)
        if (length(grep("^ERCC",rownames(norm$data))) > 0) {
                erccs = grep("^ERCC",rownames(norm$data));
                points(xes[erccs],p[erccs],col="black",cex=2)
        }
	return(list(K = krt,P = p, S = s, data = norm$data, labels = norm$labels, Kerr = Kerr, weights=norm$weights));
}
