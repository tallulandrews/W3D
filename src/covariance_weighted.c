void covariance_weighted (double* x, double* wx, double* y, double* wy, int* n, double* covar) {
	covar[0] = (double) *n;
	double sum_x = 0.0;
	double sum_wx = 0.0;
	double sum_y = 0.0;
	double sum_wy = 0.0;
	int i;
	for (i = 0; i < *n; i++) {
		sum_x = sum_x + wx[i]*x[i];
		sum_wx = sum_wx + wx[i];
		sum_y = sum_y + wy[i]*y[i];
		sum_wy = sum_wy + wy[i];
	}
	double sum_w = 0.0;
	double sum_cov = 0.0;
	for (i = 0; i < *n; i++) {
		double w = wx[i]*wy[i];
		sum_w = sum_w + w;
		sum_cov = sum_cov + w*(x[i]-sum_x/sum_wx)*(y[i]-sum_y/sum_wy);
	}
	covar[0] = 1.0/(sum_w-1.0)*sum_cov;
}
