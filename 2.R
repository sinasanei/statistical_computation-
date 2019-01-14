data=c(15.98247,16.47688,17.2935,19.29195,15.62418,19.17946,19.83762,17.66685,17.52842,14.2797,15.65064,15.46015,17.78624,16.52783,18.20742,16.9908,17.93268,21.80852,16.5106,18.25084,19.66132,15.68794,17.6264,15.05916,15.99123,16.53632,12.70026,16.52058,18.89823,16.33826,16.92829,17.40366,16.97416,15.52522,18.56712,17.70126,18.35067,14.88808,17.96294,16.6413,18.52277,17.60614,18.28284,17.21295,17.11902,18.32105,13.26564,16.90886,18.00639,17.81292,16.91042,18.82537,16.75111,15.87344,14.4237,14.79683,16.22904,17.07457,17.67224,16.62292,19.34927,16.12205,16.8359,16.30263,17.62287,15.9437,16.91413,18.1876,14.61401,18.94885,16.33257,18.65407,16.36635,16.30884,16.90533,19.1134,18.85181,16.55264,18.25012,20.19427,16.73713,17.90763,16.59451,16.27072,18.13678,15.63045,17.90084,15.0235,15.87718,15.21189,15.84529,14.22941,17.58531,18.95757,18.25935,18.36956,16.82062,16.63638,18.45567,17.42606)
hist(data, breaks = 10)
summary(data)

x=data
for ( i in 1:100){ 
  n_r ( x, 0.00001 , i)}


mle_nr=function(xvec,stop_crit){
  init=median(xvec)-1;
  n=length(xvec);
  nn = 0
  mu_seq = NULL
  mu_current=init;
  ####compute first derivative of log-likelihood #####
  first_d=n-2*sum(exp(-xvec+mu_current)/(1+exp(-xvec+mu_current)));
  ### Continue algorithm until the first derivative ###
  ### of the log-likelihood is within stop criterion ##
  while(abs(first_d)>stop_crit){
    ####compute second derivative of log-likelihood #####
    nn = nn+1
    sRate= 0.001
    mu_current=-2*sum(exp(-xvec+mu_current)/(1+exp(-xvec+mu_current))^2);
    #### Newton-Raphson’s  update of estimate of mu ####
    mu_new=mu_current- sRate * (first_d/mu_current)
    mu_seq = c(mu_seq, mu_new)
    mu_current = mu_new
    # Compute first derivative of log likelihood
    first_d=n-2*sum(exp(-xvec+mu_current)/(1+exp(-xvec+mu_current)));
  }
  return(list(mu_hat=mu_current, iterations =nn, sequence = mu_seq))
  #print(nn)
}
a=mle_nr(x, 0.0001)
a
mle_nr=function(xvec,toler){
  startvalue=median(xvec);
  n=length(xvec);
  mu_curr=startvalue;
  nn=0
  mu_seq = NULL
  ####compute first derivative of log-likelihood #####
  first_derivll=n-2*sum(exp(-xvec+mu_curr)/(1+exp(-xvec+mu_curr)));
  ### Continue algorithm until the first derivative ###
  ### of the log-likelihood is within stop criterion ##
  while(abs(first_derivll)>toler){
    ####compute second derivative of log-likelihood #####
    second_derivll=-2*sum(exp(-xvec+mu_curr)/(1+exp(-xvec+mu_curr))^2);
    #### Newton-Raphson’s  update of estimate of mu ####
    mu_new=mu_curr-first_derivll/second_derivll;
    mu_seq = c(mu_seq, mu_new);
    mu_curr=mu_new;
    ####compute first derivative of log-likelihood #####
    first_derivll=n-2*sum(exp(-xvec+mu_curr)/(1+exp(-xvec+mu_curr)));
    nn=nn+1
  }
  return (list(thetahat=mu_curr, iterations =nn ,sequence = mu_seq))
}
mle_nr(x,00.0000000001)
