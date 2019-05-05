# Profile-Bayesian

This folder contains the source code of the manuscript "Profile Bayesian estimation for extrapolationof historical adult data to pediatric drugdevelopment".

Under continuous endpoints:

The scenarios considered in the R script are:

Under null hypothesis:

H0-0.       Under H0: (mu_a=1, var_a=10^2, mu_p=-0.05, var_p=5^2)<br/>
H0-1.       Under H0: (mu_a=1, var_a=10^2, mu_p=0, var_p=5^2)<br/>
H0-2.       Under H0: (mu_a=1, var_a=10^2, mu_p=0, var_p=10^2)<br/>
H0-3.       Under H0: (mu_a=1, var_a=10^2, mu_p=0, var_p=15^2)<br/>
H0-4.       Under H0: (mu_a=1, var_a=10^2, mu_p=0, var_p=20^2)<br/>

Under alternative hypothesis:

H1-0.       Under H1: (mu_a=1, var_a=10^2, mu_p=0.5, var_p=5^2)<br/>
H1-1.       Under H1: (mu_a=1, var_a=10^2, mu_p=0.5, var_p=10^2)<br/>
H1-2.       Under H1: (mu_a=1, var_a=10^2, mu_p=0.5, var_p=15^2)<br/>
H1-3.       Under H1: (mu_a=1, var_a=10^2, mu_p=0.5, var_p=20^2)<br/>
H1-4.       Under H1: (mu_a=1, var_a=10^2, mu_p=0.8, var_p=5^2)<br/>
H1-5.       Under H1: (mu_a=1, var_a=10^2, mu_p=0.8, var_p=10^2)<br/>
H1-6.       Under H1: (mu_a=1, var_a=10^2, mu_p=0.8, var_p=15^2)<br/>
H1-7.       Under H1: (mu_a=1, var_a=10^2, mu_p=0.8, var_p=20^2)<br/>
H1-8.       Under H1: (mu_a=1, var_a=10^2, mu_p=1, var_p=10^2)<br/>
H1-9.       Under H1: (mu_a=1, var_a=10^2, mu_p=1.5, var_p=10^2)<br/>
H1-10.       Under H1: (mu_a=1, var_a=10^2, mu_p=2, var_p=10^2)<br/>

The scenarios considered in the manuscript are (will be updated):
---
#H0-3.       Under H0: (mu_a=1, var_a=10^2, mu_p=0, var_p=10^2)<br/>
#H1-1.       Under H1: (mu_a=1, var_a=10^2, mu_p=0.5, var_p=20^2)<br/>
#H1-2.       Under H1: (mu_a=1, var_a=10^2, mu_p=0.8, var_p=20^2)<br/>
#H1-3.       Under H1: (mu_a=1, var_a=10^2, mu_p=0.8, var_p=15^2)<br/>
#H1-4.       Under H1: (mu_a=1, var_a=10^2, mu_p=0.8, var_p=10^2)<br/>
#H1-5.       Under H1: (mu_a=1, var_a=10^2, mu_p=1, var_p=10^2)<br/>
#H1-6.       Under H1: (mu_a=1, var_a=10^2, mu_p=1.5, var_p=10^2)<br/>
---
In addition, illustration plot for rejection region and hypothetical power are going to be included.

###5/2/2019###

Four outputs are added:<br/>

conditional_dataout_half_19: condition: hat(mu_p)>1/2*hat(mu_a), with 0.1 on non-informative prior<br/>
conditional_dataout_half_55: condition: hat(mu_p)>1/2*hat(mu_a), with 0.5 on non-informative prior<br/>
conditional_dataout_one_19: condition: hat(mu_p)>1*hat(mu_a), with 0.1 on non-informative prior<br/>
conditional_dataout_one_55: condition: hat(mu_p)>1*hat(mu_a), with 0.5 on non-informative prior<br/>
