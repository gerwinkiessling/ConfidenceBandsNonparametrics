###Oxford Application Writing Sample


install.packages("expm")
library(expm)

install.packages("ggplot2")
library(ggplot2)

install.packages("pracma")
library(pracma)

install.packages("sfsmisc")
library(sfsmisc)

install.packages("locpol")
library(locpol)

install.packages('latex2exp')
library(latex2exp)

install.packages("readxl")
library(readxl)

### Simulations
#I begin by writing some funcitons. 


epa_kernel<- function(x){
  #The Epanechnikov kernel
  
  ##arguments: x: point where kernel is evaluated
  ##returns: K: value of the kernel function
  K<-3/4*(1-x^2)*(-1<=x)*(x<=1)
  return(K)
}

X<-function(x, eval, p){
  # regressor matrix function
  
  ##arguments: -x: vector of data points
  #            -eval: fixed x
  #            -p: degree for the polynomial
  
  ##returns:   -reg_mat: regressor matrix
  reg_mat<-matrix(NA, nrow=length(x), ncol=p+1)
  for (i in 1:(p+1)){
    reg_mat[,i] <- (x-eval)^(i-1)
  }
  return(reg_mat)
}


kernel_matrix<-function(x, eval, h){
  #diagonal W matrix with kernels
  
  ##arguments: -x: vector of data points
  #            - eval: fixed x for evaluation
  #            - h: bandwidth
  
  ##returns:   - K: kernel matrix
  K<-diag(epa_kernel((x-eval)/h))
  return(K)
}

del_hat<-function(x, x_eval, h, p){
  #function for performing local linear regression
  
  #arguments: - x: vector of data points
  #           - x_eval: fixed x for evaluation
  #           - h: bandwidth
  #           - p: degree of polynomial
  del<-solve(crossprod(t(crossprod(X(x, x_eval, p), kernel_matrix(x, x_eval, h))),X(x, x_eval, p)), tol=1e-50)%*%crossprod(X(x, x_eval, p), kernel_matrix(x, x_eval, h))
  return(del)
}

kernel_deriv<-function(x, xeval, h){
  #The kernel derivative function. Since I evaluate the derivative at a finite number of grid points, I never evaluate it
  # at the point where it is not defined. In a second version of the paper, I implement the simulation using the 
  #smooth Gaussian kernel. Results are available upon request.
  
  ##arguments: - x: data point
  #            - xeval: fixed x
  #            - h: bandwidth
  
  ##returns: - deriv: derivative
  deriv<-(3/4)*(-2*((x-xeval)/h))*(-1/h)*(abs(((x-xeval)/h))<1)
  return(deriv)
}

calc_S1_S2_func<-function(x, xeval, h){
  #calculating terms using Wassermann's (2007) S1, S2 representation
  
  #arguments: - x: vector of data points
  #           - xeval: fixed x
  #           - h: bandwidth
  
  S2_ind<-c()
  for (j in 1:length(x)){
    S2_ind[j]<-epa_kernel((x[j]-xeval)/h)*(x[j]-xeval)^2
  }
  S2<-sum(S2_ind)
  
  S1_ind<-c()
  for (j in 1:length(x)){
    S1_ind[j]<-epa_kernel((x[j]-xeval)/h)*(x[j]-xeval)
  }
  S1<-sum(S1_ind)
  
  b<-c()
  for (i in 1:length(x)){
    b[i] <- epa_kernel((x[i]-xeval)/h)*(S2-(x[i]-xeval)*S1)
  }
  
  return(b)
}


draw_errors_function<-function(x, design, n, mean_error_norm, sd_error_norm, shape_p, scale_p, hetero_1, hetero_2, deg_freedom){
  # function with different error designs
  
  ##arguments: - x: vector of data points
  #            - design: string for error design
  #            - n: sample size
  #            - mean_error_norm: mean of the normal distribution
  #            - sd_error_norm: standard error of normal distribution
  #            - shape_p: shape parameter of gamma distribution
  #            - scale_p: scale parameter of gamma distribution
  #            - hetero_1: parameter 1 for heteroskeasticity
  #            - hetero_2: parameter 2 for heteroskedasticity#
  #            - deg_freedom: degrees of freedom of t-distribution
  
  ##returns:   - errors: vector of errors
  if (design=="iidn"){
    errors<-rnorm(n, mean=mean_error_norm, sd=sd_error_norm)
  } else if (design=="skewed"){
    errors <- rgamma(n, shape=shape_p, scale=scale_p)-shape_p*scale_p
  }
  else if (design=="heteroskedastic"){
    errors<-rnorm(n, mean=mean_error_norm, sd=sd_error_norm)*sqrt(hetero_1+hetero_2*x)
  }
  else if (design=="heavy_tails"){
    errors<-rt(n, deg_freedom)
  }
  return(errors)
}



cv_loc_lin<-function(x, Y, h){
  #performs least-squares cross-validation
  
  ##arguments: - x: vector with x values
  #            - Y: vector with Y values
  #            - h: vector with bandwidths
  
  ##returns: - h_lin_opt: optimal bandwidth
  imse_lin<-c()
  opt_linear<-matrix(NA, nrow=length(x), ncol=length(h))
  for (m in 1:length(h)){
    L_ii<-c()
    gh_x<-c()
    for (k in 1:length(x)){
      xeval=x[k]
      b<-calc_S1_S2_func(x=x, xeval=xeval, h=h[m])
      l<-c()
      for (i in 1:length(x)){
        l[i]<-b[i]/sum(b)
        if(i==k){
          L_ii[k]=l[i]
        }
      }
      gh_x[k]<-l%*%Y
    }
    imse_lin[m] <- mean(((Y-gh_x)/(1-L_ii))^2)
    opt_linear[,m]<-gh_x
  }
  index<-which.min(imse_lin)
  h_lin_opt=h[index]
  return(h_lin_opt=h_lin_opt)
}

gh_x_function<-function(n, p1, p2, ncp, x_seq, design, mean_error_norm, sd_error_norm, shape_p, scale_p, 
                       hetero_1, hetero_2, undersmooth, factor_undersmooth, deg_freedom){
  #performs local linear regression
  
  ##arguments: - n: sample size, 
  #            - p1: parameter 1 of beta distribution
  #            - p2: parameter 2 of beta distribution
  #            - ncp: non-centrality parameter of beta distribution
  #            - x_seq: sequence of x_values 
  #            - design: error design 
  #            - mean_error_norm: mean of normally distributed errors
  #            - sd_error_norm: standard error of normally distributed errors 
  #            - shape_p: shape parameter of gamma distribution 
  #            - scale_p: scale parameter of gamma distribution
  #            - hetero_1: parameter 1 for heteroskedasticity
  #            - hetero_2: parameter 2 for heteroskedasticity
  #            - undersmooth: boolean for undersmoothing
  #            - factor_undersmooth: amount of undersmoothing
  #            - deg_freedom: degrees of freedom of t-distribtion
  
  ##returns: - x: x values
  #          - Y: Y values, 
  #          - gh_x: estimated g, 
  #          - gh_x_i: estimated g for xi, 
  #          - trace_L: trace of smoothing matrix 
  #          - vtilde=vtilde: trace of product of smoothing matrix
  #          - h: bandwidth, 
  #          - smoothing_matrix_grid: smoothing matrix with xi

  
  x<-rbeta(n, p1, p2, ncp=ncp)
  Y<-x*(sin(3*pi*x))^3+draw_errors_function(x, design, n, mean_error_norm, sd_error_norm, shape_p, scale_p, hetero_1, hetero_2, deg_freedom)
  h_thumb<- thumbBw(x, Y, 1, EpaK)
  red_h_seq <- seq(from=h_thumb-0.05, to=h_thumb+0.05, length.out=5)
  if (undersmooth=="no"){
  h<-cv_loc_lin(x, Y, red_h_seq)
  }
  else if (undersmooth=="yes"){
    h<-n^(factor_undersmooth)*cv_loc_lin(x, Y, red_h_seq)
  }
  gh_x<-c()
  gh_x_i<-c()
  smoothing_matrix <- matrix(NA, nrow=length(x), ncol=length(x))
  smoothing_matrix_grid <- matrix(NA, nrow=length(x_seq), ncol=length(x))
  for (k in 1:length(x_seq)){
    xeval=x_seq[k]
    b<-calc_S1_S2_func(x=x, xeval=xeval, h=h)
    l<-c()
    for (i in 1:length(x)){
      l[i]<-b[i]/sum(b)
    }
    smoothing_matrix_grid[k,] <- l
    gh_x[k]<-l%*%Y
  }
  for (k in 1:length(x)){
    xeval=x[k]
    b<-calc_S1_S2_func(x=x, xeval=xeval, h=h)
    l<-c()
    for (i in 1:length(x)){
      l[i]<-b[i]/sum(b)
    }
    smoothing_matrix[k,] <- l
    gh_x_i[k]<-l%*%Y
  }
  trace_L<-sum(diag(smoothing_matrix))
  vtilde<-sum(diag(crossprod(smoothing_matrix, smoothing_matrix)))
  return(list(x=x, Y=Y, gh_x=gh_x, gh_x_i=gh_x_i, trace_L=trace_L, vtilde=vtilde, h=h, smoothing_matrix_grid=smoothing_matrix_grid))
}


error_function<-function(x, Y, gh_x_i, smoothing_matrix_grid){
  #calculating standard error of g
  
  #arguments: - x: x values
  #           - Y: y values
  #           - gh_x_i: estimated values of regression function at xi
  #           - smoothing_matrix_grid: smoothing matrix at xi values
  #returns:   - s_x: standard error of g hat
  
  Z<- (Y-gh_x_i)^2
  h<- thumbBw(x, Z, 1, EpaK)
    sigma_x_i<-c()
     for (k in 1:length(x)){
     xeval=x[k]
     b<-calc_S1_S2_func(x=x, xeval=xeval, h=h)
     l<-c()
     for (i in 1:length(x)){
       l[i]<-b[i]/sum(b)
     }
     sigma_x_i[k]<-l%*%Z
     }
    s_x<-as.vector(sigma_x_i%*%t(smoothing_matrix_grid^2))
    s_x[s_x<0] <- 0
    s_x <- sqrt(s_x)
   
  return(s_x)
}

derivative_function<-function(xeval, h, x){
  # calculates integrand to get kappa
  
  # arguments: - xeval: fixed x
  #            - h: bandwidth
  #            - x: x values
  # returns: - T_norm: integrand
  #           - l_norm: norm of smoother
  S2_ind<-c()
  for (j in 1:length(x)){
    S2_ind[j]<-epa_kernel((x[j]-xeval)/h)*(x[j]-xeval)^2
  }
  S2<-sum(S2_ind)
  
  S1_ind<-c()
  for (j in 1:length(x)){
    S1_ind[j]<-epa_kernel((x[j]-xeval)/h)*(x[j]-xeval)
  }
  S1<-sum(S1_ind)
  #taking care of bi
  b<-c()
  for (i in 1:length(x)){
    b[i] <- epa_kernel((x[i]-xeval)/h)*(S2-(x[i]-xeval)*S1)
  }
  l<-c()
  for (i in 1:length(x)){
    l[i]<-b[i]/sum(b)
  }
  l_norm<-sqrt(sum(l^2))
  
  
  first_part<-c()
  for (i in 1:length(x)){
    first_part[i]<-kernel_deriv(x[i], xeval, h)*(S2-(x[i]-xeval)*S1)
  }
  
  derivS2_part1_ind<-c()
  for (j in 1:length(x)){
    derivS2_part1_ind[j]<-kernel_deriv(x[j], xeval, h)*(x[j]-xeval)^2
  }
  derivS2_part1<-sum(derivS2_part1_ind)
  derivS2_part2_ind<-c()
  for (j in 1:length(x)){
    derivS2_part2_ind[j]<-epa_kernel((x[j]-xeval)/h)*(-2*(x[j]-xeval))
  }
  derivS2_part2<-sum(derivS2_part2_ind)
  
  derivS1_part1_ind<-c()
  for (j in 1:length(x)){
    derivS1_part1_ind[j]<-epa_kernel(((x[j]-xeval)/h))*(x[j]-xeval)
  }
  derivS1_part1<-sum(derivS1_part1_ind)
  
  derivS1_part2_ind<-c()
  for (j in 1:length(x)){
    derivS1_part2_ind[j]<-kernel_deriv(x[j], xeval, h)*(x[j]-xeval)-epa_kernel((x[j]-xeval)/h)
  }
  derivS1_part2<-sum(derivS1_part2_ind)
  
  second_part<-c()
  for (i in 1:length(x)){
    second_part[i] <- epa_kernel((x[i]-xeval)/h)*(derivS2_part1+derivS2_part2+derivS1_part1-(x[i]-xeval)*derivS1_part2)
  }
  del_b<-c()
  for (i in 1:length(x)){
    del_b[i]<-first_part[i]+second_part[i]
  }
  
  del_T<-c()
  for (i in 1:length(x)){
    del_T[i] <- (del_b[i]*sqrt(sum(b^2))-b[i]*(2*sum(del_b)/(sqrt(sum(b^2)))))/sum(b^2)
  }
  T_norm<-sqrt(sum(del_T^2))
  return(list(T_norm=T_norm, l_norm=l_norm))
}

f<-function(n, alpha, seq, kappa, trace_L){
  # calculates c
  ##arguments: n: sample size
  #            alpha: critical values 
  #            seq: sequence of c values 
  #            kappa: kappa0 value 
  #            trace_L: trace of the smoothing matrix
  ##returns: solving_c: c that solves the equation
  #          difference: difference between LHS and RHS
  m<-n-trace_L
  
  difference<-c()
  difference<-abs((1-pt(q=seq, df=m))+(kappa/pi)*(1+(seq)^2/m)^(-m/2)-alpha)
  index=which.min(difference)
  solving_c=seq[index]
  return(list(solving_c=solving_c, difference=difference))
}


uniform_ll<-function(n, p1, p2, ncp, x_seq, lower_bound_c, upper_bound_c, steps_c, alpha, design, 
                     mean_error_norm, sd_error_norm, shape_p, scale_p, 
                     heteroskedasticity_robust, hetero_1, hetero_2, undersmooth, factor_undersmooth, 
                     deg_freedom){
  
  #calculates uniform confidence bands
  
  ##arguments: - n: sample size, 
  #            - p1: parameter 1 of beta distribution
  #            - p2: parameter 2 of beta distribution
  #            - ncp: non-centrality parameter of beta distribution
  #            - x_seq: sequence of x_values 
  #            - lower_bound_c: lower bound of c sequence
  #            - upper bound of c sequence
  #            - steps: steps at which to evaluate c
  #            - alpha: critical value
  #            - design: error design 
  #            - mean_error_norm: mean of normally distributed errors
  #            - sd_error_norm: standard error of normally distributed errors 
  #            - shape_p: shape parameter of gamma distribution 
  #            - scale_p: scale parameter of gamma distribution
  #            - hetero_1: parameter 1 for heteroskedasticity
  #            - hetero_2: parameter 2 for heteroskedasticity
  #            - undersmooth: boolean for undersmoothing
  #            - factor_undersmooth: amount of undersmoothing
  #            - deg_freedom: degrees of freedom of t-distribtion
  
  ##returns:   - all_in: boolean whether function is covered 
  #            - c: c for confidence band 
  #            - results_df: data frame with results
  #            - x: x values
  #            - Y: Y values
  #            - h: bandwidth
  
  
  g_x<-x_seq*(sin(3*pi*x_seq))^3
  
  results<-gh_x_function(n, p1, p2, ncp, x_seq, design, mean_error_norm, sd_error_norm, shape_p, scale_p, hetero_1, hetero_2, 
                        undersmooth, factor_undersmooth, deg_freedom)
  
  results_variance<-data.frame('x'=results$x, 'gh_x_i'=results$gh_x_i, 'Y'=results$Y)
  results_df<-data.frame('gh_x'=results$gh_x, 'g_x'=g_x)
  trace_L<- results$trace_L
  
  l_norm<-c()
  points<-c()
  
  for (i in 1:length(x_seq)){
    results_derivative <- derivative_function(x_seq[i], results$h, results$x)
    l_norm[i]<-results_derivative$l_norm
    points[i]<-results_derivative$T_norm
  }
  
  kappa0<-integrate.xy(x_seq, points, min(x_seq), max(x_seq))
  
  c_seq=seq(from=lower_bound_c, to=upper_bound_c, by=steps_c)
  
  result_c<-f(alpha=alpha, seq=c_seq, kappa=kappa0, trace_L=trace_L , n=n)
  c<-result_c$solving_c
  
  if (heteroskedasticity_robust=="no"){
  s<-sqrt(sum((results$Y-results$gh_x_i)^2)/(n-2*trace_L+results$vtilde))*l_norm
  }
  
  else if (heteroskedasticity_robust=="yes"){
  s<-error_function(results$x, results$Y, results$gh_x_i, results$smoothing_matrix_grid)
  }
  
  
  CI_lower<-results_df$gh_x-c*s
  CI_upper<-results_df$gh_x+c*s
  if(length(CI_lower)>0){
  results_df<-cbind(results_df, CI_lower, CI_upper)

  
  
  boolean<-(results_df$CI_lower<results_df$g_x&results_df$g_x<results_df$CI_upper)
  
  all_in<-all(boolean)
}
else{
  all_in<-NA
  results_df<-data.frame()
}
  
  return(list(all_in=all_in, c=c, results_df=results_df, x=results$x, Y=results$Y, h=results$h))
}

replication_function<-function(reps, n, p1, p2, ncp, lower_bound_x, upper_bound_x, steps_x, 
                               lower_bound_c, upper_bound_c, steps_c, alpha, design,
                               mean_error_norm, sd_error_norm, shape_p, scale_p, 
                               heteroskedasticity_robust, hetero_1, hetero_2, undersmooth, factor_undersmooth, 
                               deg_freedom){
  
  #performs uniform_ll function reps times
  
  ##returns: plot: ggplot
  #          uniform_coverage: coverage rate
  #          solving c: c used in confidence band
  #          solving_h: optimal bandwidth

  

  
  x_seq=seq(from=lower_bound_x, upper_bound_x, steps_x)
  
  
  replication<-replicate(reps, uniform_ll(n=n, p1=p1, p2=p2, ncp=ncp, x_seq,  
                                          lower_bound_c=lower_bound_c, upper_bound_c=upper_bound_c, steps_c=steps_c, alpha=alpha, 
                                          design=design, mean_error_norm=mean_error_norm, sd_error_norm=sd_error_norm, 
                                          shape_p=shape_p, scale_p=scale_p, heteroskedasticity_robust, hetero_1, hetero_2, 
                                          undersmooth, factor_undersmooth, deg_freedom))
  uniform_coverage<-mean(as.numeric(replication['all_in',]), na.rm=TRUE)
  solving_c<-as.numeric(replication['c',])
  solving_h<-as.numeric(replication['h',])
  
  
  df<-cbind(x_seq, data.frame(replication['results_df',1]))
  df_points<-data.frame('x'=replication['x',reps], 'Y'=replication['Y',1])
  
  plot<-ggplot(df, aes(x=x_seq)) +
    geom_line(aes(y=results_df.g_x, colour="True regression function")) +
    geom_line(aes(y=results_df.gh_x, colour="Estimated regression function")) +
    geom_ribbon(aes(ymin=results_df.CI_lower,ymax=results_df.CI_upper, fill="Confidence band"),alpha=0.3)+
    geom_point(data=df_points, aes(x=x,y=Y), color="black", size=0.05) + 
    theme_bw() + xlab(TeX('$x$')) + 
    scale_fill_manual("",values="grey12") +
    scale_colour_manual("", 
                        breaks = c("True regression function", "Estimated regression function"),
                        values = c("red", "black")) +
    theme(legend.position="bottom", legend.margin=margin(-10, 0, 0, 0), legend.direction="vertical", 
          axis.title.y = element_blank(), text = element_text(size = 15))
  
  return(list(plot=plot, uniform_coverage=uniform_coverage, solving_c=solving_c, solving_h=solving_h))
}

set.seed(1)

sim1<-replication_function(reps=1000,  n=200, p1=1.1, p2=1.1, ncp=0, lower_bound_x=0.01, upper_bound_x = 0.01, steps_x=0.01,
                     lower_bound_c=2, upper_bound_c = 7, steps_c = 0.01, alpha=0.01, 
                     design="iidn", mean_error_norm = 0, sd_error_norm=0.5, shape_p=1.2, scale_p=0.45, 
                     heteroskedasticity_robust = "no", hetero_1=0.1, hetero_2=1, undersmooth="no", factor_undersmooth = -1/20, 
                     deg_freedom = 4)
sim1
set.seed(2)

sim2<-replication_function(reps=1000,  n=200, p1=1.1, p2=1.1, ncp=0, lower_bound_x=0.01, upper_bound_x = 0.09, steps_x=0.01,
                           lower_bound_c=2, upper_bound_c = 6, steps_c = 0.01, alpha=0.01, 
                           design="skewed", mean_error_norm = 0, sd_error_norm=0.5, shape_p=1.2, scale_p=0.45,
                           heteroskedasticity_robust = "no", hetero_1=0.1, hetero_2=1, undersmooth="no", factor_undersmooth = -1/20, 
                           deg_freedom = 4)
sim2
set.seed(3)

sim3<-replication_function(reps=1000,  n=200, p1=1.1, p2=1.1, ncp=0, lower_bound_x=0.01, upper_bound_x = 0.99, steps_x=0.01,
       lower_bound_c=2, upper_bound_c = 6, steps_c = 0.01, alpha=0.01, 
       design="heteroskedastic", mean_error_norm = 0, sd_error_norm=0.5, shape_p=1.2, scale_p=0.45,
       heteroskedasticity_robust = "yes", hetero_1=0.1, hetero_2=1, undersmooth="no", factor_undersmooth = -1/20, 
       deg_freedom = 4)
sim3
set.seed(4)

sim4<-replication_function(reps=1000,  n=200, p1=1.1, p2=1.1, ncp=0, lower_bound_x=0.01, upper_bound_x = 0.99, steps_x=0.01,
                           lower_bound_c=2, upper_bound_c = 6, steps_c = 0.01, alpha=0.01, 
                           design="heteroskedastic", mean_error_norm = 0, sd_error_norm=0.5, shape_p=1.2, scale_p=0.45,
                           heteroskedasticity_robust = "no", hetero_1=0.1, hetero_2=1, undersmooth="no", factor_undersmooth = -1/20, 
                           deg_freedom = 4)
sim4
set.seed(5)

sim5<-replication_function(reps=1000,  n=200, p1=1.1, p2=1.1, ncp=0, lower_bound_x=0.01, upper_bound_x = 0.99, steps_x=0.01,
                           lower_bound_c=2, upper_bound_c = 6, steps_c = 0.01, alpha=0.01, 
                           design="iidn", mean_error_norm = 0, sd_error_norm=0.5, shape_p=1.2, scale_p=0.45, 
                           heteroskedasticity_robust = "no", hetero_1=0.1, hetero_2=1, undersmooth="yes", factor_undersmooth = -1/20, 
                           deg_freedom = 4)
sim5
set.seed(6)

sim6<-replication_function(reps=1000,  n=200, p1=1.1, p2=1.1, ncp=0, lower_bound_x=0.01, upper_bound_x = 0.99, steps_x=0.01,
                           lower_bound_c=2, upper_bound_c = 6, steps_c = 0.01, alpha=0.01, 
                           design="skewed", mean_error_norm = 0, sd_error_norm=0.5, shape_p=1.2, scale_p=0.45, 
                           heteroskedasticity_robust = "no", hetero_1=0.1, hetero_2=1, undersmooth="yes", factor_undersmooth = -1/20, 
                           deg_freedom = 4)
sim6
set.seed(77)

sim7<-replication_function(reps=1000,  n=200, p1=1.1, p2=1.1, ncp=0, lower_bound_x=0.01, upper_bound_x = 0.99, steps_x=0.01,
                           lower_bound_c=2, upper_bound_c = 6, steps_c = 0.01, alpha=0.01, 
                           design="heteroskedastic", mean_error_norm = 0, sd_error_norm=0.5, shape_p=1.2, scale_p=0.45,
                           heteroskedasticity_robust = "yes", hetero_1=0.1, hetero_2=1, undersmooth="yes", factor_undersmooth = -1/20, 
                           deg_freedom = 4)
sim7
set.seed(8)

sim8<-replication_function(reps=1000,  n=200, p1=1.1, p2=1.1, ncp=0, lower_bound_x=0.01, upper_bound_x = 0.99, steps_x=0.01,
                           lower_bound_c=2, upper_bound_c = 6, steps_c = 0.01, alpha=0.01, 
                           design="heteroskedastic", mean_error_norm = 0, sd_error_norm=0.5, shape_p=1.2, scale_p=0.45, 
                           heteroskedasticity_robust = "no", hetero_1=0.1, hetero_2=1, undersmooth="yes", factor_undersmooth = -1/20, 
                           deg_freedom = 4)
sim8
set.seed(9)

sim9<-replication_function(reps=1000, n=400, p1=1.1, p2=1.1, ncp=0, lower_bound_x=0.01, upper_bound_x = 0.99, steps_x=0.01,
                           lower_bound_c=2, upper_bound_c = 6, steps_c = 0.01, alpha=0.01, 
                           design="iidn", mean_error_norm = 0, sd_error_norm=0.5, shape_p=1.2, scale_p=0.45, 
                           heteroskedasticity_robust = "no", hetero_1=0.1, hetero_2=1, undersmooth="no", factor_undersmooth = -1/20, 
                           deg_freedom = 4)
sim9$plot
ggsave(filename="sim9plot.png", plot=sim9$plot, device="png", dpi=500)

set.seed(10)

sim10<-replication_function(reps=1000, n=400, p1=1.1, p2=1.1, ncp=0, lower_bound_x=0.01, upper_bound_x = 0.99, steps_x=0.01,
                           lower_bound_c=2, upper_bound_c = 6, steps_c = 0.01, alpha=0.01, 
                           design="skewed", mean_error_norm = 0, sd_error_norm=0.5, shape_p=1.2, scale_p=0.45,
                           heteroskedasticity_robust = "no", hetero_1=0.1, hetero_2=1, undersmooth="no", factor_undersmooth = -1/20, 
                           deg_freedom = 4)

set.seed(11)

sim11<-replication_function(reps=1000, n=400, p1=1.1, p2=1.1, ncp=0, lower_bound_x=0.01, upper_bound_x = 0.99, steps_x=0.01,
                           lower_bound_c=2, upper_bound_c = 6, steps_c = 0.01, alpha=0.01, 
                           design="heteroskedastic", mean_error_norm = 0, sd_error_norm=0.5, shape_p=1.2, scale_p=0.45, 
                           heteroskedasticity_robust = "yes", hetero_1=0.1, hetero_2=1, undersmooth="no", factor_undersmooth = -1/20, 
                           deg_freedom = 4)

ggsave(filename="sim11plot.png", plot=sim11$plot, device="png", dpi=500)

dev.off()
set.seed(12)

sim12<-replication_function(reps=1000, n=400, p1=1.1, p2=1.1, ncp=0, lower_bound_x=0.01, upper_bound_x = 0.99, steps_x=0.01,
                           lower_bound_c=2, upper_bound_c = 6, steps_c = 0.01, alpha=0.01, 
                           design="heteroskedastic", mean_error_norm = 0, sd_error_norm=0.5, shape_p=1.2, scale_p=0.45, 
                           heteroskedasticity_robust = "no", hetero_1=0.1, hetero_2=1, undersmooth="no", factor_undersmooth = -1/20, 
                           deg_freedom = 4)
sim12

set.seed(13)

sim13<-replication_function(reps=1000, n=400, p1=1.1, p2=1.1, ncp=0, lower_bound_x=0.01, upper_bound_x = 0.99, steps_x=0.01,
                           lower_bound_c=2, upper_bound_c = 6, steps_c = 0.01, alpha=0.01, 
                           design="iidn", mean_error_norm = 0, sd_error_norm=0.5, shape_p=1.2, scale_p=0.45, 
                           heteroskedasticity_robust = "no", hetero_1=0.1, hetero_2=1, undersmooth="yes", factor_undersmooth = -1/20, 
                           deg_freedom = 4)

set.seed(14)

sim14<-replication_function(reps=1000, n=400, p1=1.1, p2=1.1, ncp=0, lower_bound_x=0.01, upper_bound_x = 0.99, steps_x=0.01,
                           lower_bound_c=2, upper_bound_c = 6, steps_c = 0.01, alpha=0.01, 
                           design="skewed", mean_error_norm = 0, sd_error_norm=0.5, shape_p=1.2, scale_p=0.45, 
                           heteroskedasticity_robust = "no", hetero_1=0.1, hetero_2=1, undersmooth="yes", factor_undersmooth = -1/20, 
                           deg_freedom = 4)

set.seed(15)

sim15<-replication_function(reps=1000, n=400, p1=1.1, p2=1.1, ncp=0, lower_bound_x=0.01, upper_bound_x = 0.99, steps_x=0.01,
                           lower_bound_c=2, upper_bound_c = 6, steps_c = 0.01, alpha=0.01, 
                           design="heteroskedastic", mean_error_norm = 0, sd_error_norm=0.5, shape_p=1.2, scale_p=0.45, 
                           heteroskedasticity_robust = "yes", hetero_1=0.1, hetero_2=1, undersmooth="yes", factor_undersmooth = -1/20, 
                           deg_freedom = 4)

sim15
set.seed(16)

sim16<-replication_function(reps=1000, n=400, p1=1.1, p2=1.1, ncp=0, lower_bound_x=0.01, upper_bound_x = 0.99, steps_x=0.01,
                           lower_bound_c=2, upper_bound_c = 6, steps_c = 0.01, alpha=0.01, 
                           design="heteroskedastic", mean_error_norm = 0, sd_error_norm=0.5, shape_p=1.2, scale_p=0.45, 
                           heteroskedasticity_robust = "no", hetero_1=0.1, hetero_2=1, undersmooth="yes", factor_undersmooth = -1/20, 
                           deg_freedom = 4)
sim16

set.seed(17)

sim17<-replication_function(reps=1000, n=400, p1=1.1, p2=1.1, ncp=0, lower_bound_x=0.01, upper_bound_x = 0.99, steps_x=0.01,
                            lower_bound_c=2, upper_bound_c = 6, steps_c = 0.01, alpha=0.01, 
                            design="heavy_tails", mean_error_norm = 0, sd_error_norm=0.5, shape_p=1.2, scale_p=0.45, 
                            heteroskedasticity_robust = "no", hetero_1=0.1, hetero_2=1, undersmooth="no", factor_undersmooth = -1/20, 
                            deg_freedom=4)
sim17
set.seed(18)

sim18<-replication_function(reps=1000, n=400, p1=1.1, p2=1.1, ncp=0, lower_bound_x=0.01, upper_bound_x = 0.99, steps_x=0.01,
                            lower_bound_c=2, upper_bound_c = 6, steps_c = 0.01, alpha=0.01, 
                            design="heavy_tails", mean_error_norm = 0, sd_error_norm=0.5, shape_p=1.2, scale_p=0.45, 
                            heteroskedasticity_robust = "no", hetero_1=0.1, hetero_2=1, undersmooth="yes", factor_undersmooth = -1/20, 
                            deg_freedom=4)
sim18
set.seed(19)

sim19<-replication_function(reps=1000, n=200, p1=1.1, p2=1.1, ncp=0, lower_bound_x=0.01, upper_bound_x = 0.99, steps_x=0.01,
                            lower_bound_c=2, upper_bound_c = 6, steps_c = 0.01, alpha=0.01, 
                            design="heavy_tails", mean_error_norm = 0, sd_error_norm=0.5, shape_p=1.2, scale_p=0.45, 
                            heteroskedasticity_robust = "no", hetero_1=0.1, hetero_2=1, undersmooth="no", factor_undersmooth = -1/20, 
                            deg_freedom=4)
sim19
set.seed(20)

sim20<-replication_function(reps=1000, n=200, p1=1.1, p2=1.1, ncp=0, lower_bound_x=0.01, upper_bound_x = 0.99, steps_x=0.01,
                            lower_bound_c=2, upper_bound_c = 6, steps_c = 0.01, alpha=0.01, 
                            design="heavy_tails", mean_error_norm = 0, sd_error_norm=0.5, shape_p=1.2, scale_p=0.45, 
                            heteroskedasticity_robust = "no", hetero_1=0.1, hetero_2=1, undersmooth="yes", factor_undersmooth = -1/20, 
                            deg_freedom=4)
sim20


###Empirical Application

###Additional functions

empirical_gh_x_function<-function(x, Y, h_seq, x_seq){
  #empirical local linear
  
  #arguments: x: x values
  #           Y: Y values
  #           h_seq: bandwidth sequence
  #           x_seq: x_sequence
  
  #returns:  - gh_x: estimated g, 
  #          - gh_x_i: estimated g for xi, 
  #          - trace_L: trace of smoothing matrix, 
  #          - vtilde: trace of product of smoothing matrix
  #          - h: bandwidth, 
  #          - smoothing_matrix_grid: smoothing matrix with xi
  
  h<-cv_loc_lin(x, Y, h_seq)
  gh_x<-c()
  gh_x_i<-c()
  smoothing_matrix <- matrix(NA, nrow=length(x), ncol=length(x))
  smoothing_matrix_grid <- matrix(NA, nrow=length(x_seq), ncol=length(x))
  for (k in 1:length(x_seq)){
    xeval=x_seq[k]
    b<-calc_S1_S2_func(x=x, xeval=xeval, h=h)
    l<-c()
    for (i in 1:length(x)){
      l[i]<-b[i]/sum(b)
    }
    smoothing_matrix_grid[k,] <- l
    gh_x[k]<-l%*%Y
  }
  for (k in 1:length(x)){
    xeval=x[k]
    b<-calc_S1_S2_func(x=x, xeval=xeval, h=h)
    l<-c()
    for (i in 1:length(x)){
      l[i]<-b[i]/sum(b)
    }
    smoothing_matrix[k,] <- l
    gh_x_i[k]<-l%*%Y
  }
  trace_L<-sum(diag(smoothing_matrix))
  vtilde<-sum(diag(crossprod(smoothing_matrix, smoothing_matrix)))
  return(list(gh_x=gh_x, gh_x_i=gh_x_i, trace_L=trace_L, vtilde=vtilde, h=h, smoothing_matrix_grid=smoothing_matrix_grid))
}

###

plot_empirical_function<-function(x, Y, h_seq, x_seq, c_seq){
  #produces plot
  
  #arguments: x: x values
  #           Y: Y values
  #           h_seq: bandwidth sequence
  #           c_seq: sequence of c a values
  #returns:   plot: gg_plot
  
  results_empirical<-empirical_gh_x_function(x=x, Y=Y, h_seq=h_seq, x_seq=x_seq)
  s<-error_function(x, Y, results_empirical$gh_x_i, results_empirical$smoothing_matrix_grid)
  
  l_norm<-c()
  points<-c()
  
  for (i in 1:length(x_seq)){
    results_derivative <- derivative_function(x_seq[i], results_empirical$h, x)
    l_norm[i]<-results_derivative$l_norm
    points[i]<-results_derivative$T_norm
  }
  
  kappa0<-integrate.xy(x_seq, points, min(x_seq), max(x_seq))
  
  result_c<-f(alpha=alpha, seq=c_seq, kappa=kappa0, trace_L=results_empirical$trace_L , n=n)
  c<-result_c$solving_c
  
  CI_lower<-results_empirical$gh_x-c*s
  CI_upper<-results_empirical$gh_x+c*s
  
  results_df<-data.frame('x_seq'=x_seq, 'gh_x'=results_empirical$gh_x, 'CI_lower'=CI_lower, "CI_upper"=CI_upper)
  
  results_plot<-data.frame('x'=x, 'Y'=Y)
  
  
  plot<-ggplot(results_df, aes(x=x_seq)) +
    geom_line(aes(y=gh_x), color="black") +
    geom_ribbon(aes(ymin=CI_lower,ymax=CI_upper, fill="Confidence band"),alpha=0.3)+
    geom_point(data=results_plot, aes(x=x,y=Y), color="black", size=0.05) + 
    theme_bw() + xlab("Work experience (years)") + 
    scale_fill_manual("",values="grey12")+
    theme( text = element_text(size = 17), legend.position = "none")+
    ylab("Log of yearly labour income in Pesos")
  return(plot)
}

###some data cleaning

data<-read_excel("urban.xls", col_names = FALSE)
df<-data.frame(cbind(data[,6],
               data[,23],
               data[,19],
               data[,24],
               data[,25], 
               data[,16],
               data[,17], 
               data[,18]))
names <- c("wage", "hhsize", "married", 
           "bothwork", "notemployed", 
           "college", "age", 
           "female")
colnames(df)<-names
df<-df[(df[,'female']==0),]
df<-df[df[,'bothwork']==0, ]
df<-df[df[,'college']==1, ]
df<-df[df[,'wage']>0, ]
df<-df[df[,'notemployed']==0, ]
df<-df[(22<df[,'age'])&(df[,'age']<65),]
df['experience']=df['age']-22

x=as.numeric(unlist(df['experience']))
Y=log(as.numeric(unlist(df['wage'])))
index_Y<-which.max(Y)
x=x[-index_Y]
Y=Y[-index_Y]


x_seq=seq(from=min(x), to=max(x), by=1)
h_seq=seq(from=2, to=20, by=1)
c_seq<-seq(from=2, to=6, by=0.01)

plot_empirical<-plot_empirical_function(x=x, Y=Y, h_seq=h_seq, x_seq=x_seq, c_seq=c_seq)
plot_empirical
ggsave(filename="empirical_plot.png", plot=plot_empirical, device="png", dpi=500)
empirical_gh_x_function(x=x, Y=Y, h_seq=h_seq, x_seq=x_seq)
