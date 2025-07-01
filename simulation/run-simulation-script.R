# source('~/dl-survival-miRNA/simulation/aug-batch-survival-pipeline.R')
# Run Simulation ---------------------------------

################################ BE00Asso00 #####################################
## ########################### BE00Asso00_normNone ################################
# Linear risk with moderate effects
linear.moderate.out = sim.survdata(surv.data,
                                   N=10000, N.test=10000, test.ids=test.ids, sim.dataType='linear-moderate', 
                                   nonzero.genes=nonzero.genes, n_batch=10, h0=2.75, 
                                   he_train=0, he_test=0,beta_sort_train=0,beta_sort_test=0,norm_type=0,
                                   save_gene_data=T, save_surv_data=T, save_beta=T, plot_km=T, run_analysis=F,
                                   seed=1234)
linear.moderate.out[[3]] ## Look at the distribution of simulated data

#### NOTE: set save_gene_data = F from here
# Linear risk with weak effects
linear.weak.out = sim.survdata(surv.data,
                               N=10000, N.test=10000, test.ids=test.ids, sim.dataType='linear-weak', 
                               nonzero.genes=nonzero.genes, n_batch=10, h0=50, 
                               he_train=0,he_test=0,beta_sort_train=0,beta_sort_test=0,norm_type=0,
                               save_gene_data=T, save_surv_data=T, save_beta=T, plot_km=T, run_analysis=F,
                               seed=1234)
linear.weak.out[[3]]


# Squared terms
nl.qua.out = sim.survdata(surv.data,
                          N=10000, N.test=10000, test.ids=test.ids, sim.dataType='nl-quadratic', 
                          nonzero.genes=nonzero.genes, n_batch=10, h0=4.5, 
                          he_train=0,he_test=0,beta_sort_train=0,beta_sort_test=0,norm_type=0,
                          save_gene_data=T, save_surv_data=T, save_beta=T, plot_km=F, run_analysis=F,
                          seed=1234)
nl.qua.out[[3]]

# Squared terms with intercept
nl.shiftqua.out = sim.survdata(surv.data,
                          N=10000,N.test=10000,test.ids=test.ids,sim.dataType='nl-shiftquad',
                          nonzero.genes=nonzero.genes,n_batch=10,h0=21,
                          he_train=0,he_test=0,beta_sort_train=0,beta_sort_test=0,norm_type=0,
                          save_gene_data=T,save_surv_data=T,save_beta=T,plot_km=T,run_analysis=T,
                          seed=1234)
nl.shiftqua.out[[3]]

# Gene-gene interaction terms
nl.interact.out = sim.survdata(surv.data,
                               N=10000, N.test=10000, test.ids=test.ids, sim.dataType='nl-interaction', 
                               nonzero.genes=nonzero.genes, n_batch=10, h0=32, 
                               he_train=0,he_test=0,beta_sort_train=0,beta_sort_test=0,norm_type=0,
                               save_gene_data=T, save_surv_data=T, save_beta=T, plot_km=F, run_analysis=T,
                               seed=1234)
nl.interact.out[[3]]


# Sine interaction
nl.sine.out = sim.survdata(surv.data,
                           N=10000, N.test=10000, test.ids=test.ids, sim.dataType='nl-sine', 
                           nonzero.genes=nonzero.genes, n_batch=10, h0=24, 
                           he_train=0,he_test=0,beta_sort_train=0,beta_sort_test=0,norm_type=0,
                           save_gene_data=T, save_surv_data=T, save_beta=T, plot_km=F, run_analysis=T,
                           seed=1234)
nl.sine.out[[3]]




## ########################### BE00Asso00_normTC ###############################
# Linear risk with moderate effects
linear.moderate.out = sim.survdata(surv.data,
                                   N=10000, N.test=10000, test.ids=test.ids, sim.dataType='linear-moderate', 
                                   nonzero.genes=nonzero.genes, n_batch=10, h0=13, 
                                   he_train=0,he_test=0,beta_sort_train=0,beta_sort_test=0, norm_type=1,
                                   save_gene_data=T, save_surv_data=T, save_beta=T, plot_km=F, run_analysis=T,
                                   seed=1234)
linear.moderate.out[[3]] ## Look at the distribution of simulated data

#### NOTE: set save_gene_data = F from here
# Linear risk with weak effects
linear.weak.out = sim.survdata(surv.data,
                               N=10000, N.test=10000, test.ids=test.ids, sim.dataType='linear-weak', 
                               nonzero.genes=nonzero.genes, n_batch=10, h0=26, 
                               he_train=0,he_test=0,beta_sort_train=0,beta_sort_test=0,norm_type=1,
                               save_gene_data=T, save_surv_data=T, save_beta=T, plot_km=F, run_analysis=T,
                               seed=1234)
linear.weak.out[[3]]


# Squared terms
nl.qua.out = sim.survdata(surv.data,
                          N=10000, N.test=10000, test.ids=test.ids, sim.dataType='nl-quadratic', 
                          nonzero.genes=nonzero.genes, n_batch=10, h0=18, 
                          he_train=0,he_test=0,beta_sort_train=0,beta_sort_test=0,norm_type=1,
                          save_gene_data=T, save_surv_data=T, save_beta=T, plot_km=F, run_analysis=T,
                          seed=1234)
nl.qua.out[[3]]

# Squared terms with intercept
nl.shiftqua.out = sim.survdata(surv.data,
                               N=10000,N.test=10000,test.ids=test.ids,sim.dataType='nl-shiftquad',
                               nonzero.genes=nonzero.genes,n_batch=10,h0=21,
                               he_train=0,he_test=0,beta_sort_train=0,beta_sort_test=0,norm_type=1,
                               save_gene_data=T,save_surv_data=T,save_beta=T,plot_km=T,run_analysis=T,
                               seed=1234)
nl.shiftqua.out[[3]]
# Gene-gene interaction terms
nl.interact.out = sim.survdata(surv.data,
                               N=10000, N.test=10000, test.ids=test.ids, sim.dataType='nl-interaction', 
                               nonzero.genes=nonzero.genes, n_batch=10, h0=32, 
                               he_train=0,he_test=0,beta_sort_train=0,beta_sort_test=0,norm_type=1,
                               save_gene_data=T, save_surv_data=T, save_beta=T, plot_km=F, run_analysis=T,
                               seed=1234)
nl.interact.out[[3]]


# Sine interaction
nl.sine.out = sim.survdata(surv.data,
                           N=10000, N.test=10000, test.ids=test.ids, sim.dataType='nl-sine', 
                           nonzero.genes=nonzero.genes, n_batch=10, h0=24, 
                           he_train=0,he_test=0,beta_sort_train=0,beta_sort_test=0,norm_type=1,
                           save_gene_data=T, save_surv_data=T, save_beta=T, plot_km=F, run_analysis=T,
                           seed=1234)
nl.sine.out[[3]]



## ########################### BE00Asso00_normUQ ###############################
# Linear risk with moderate effects
linear.moderate.out = sim.survdata(surv.data,
                                   N=10000, N.test=10000, test.ids=test.ids, sim.dataType='linear-moderate',
                                   nonzero.genes=nonzero.genes, n_batch=10, h0=9,
                                   he_train=0,he_test=0,beta_sort_train=0,beta_sort_test=0, norm_type=2,
                                   save_gene_data=T, save_surv_data=T, save_beta=T, plot_km=F, run_analysis=T,
                                   1234)
linear.moderate.out[[3]] ## Look at the distribution of simulated data

#### NOTE: set save_gene_data = F from here
# Linear risk with weak effects
linear.weak.out = sim.survdata(surv.data,
                               N=10000, N.test=10000, test.ids=test.ids, sim.dataType='linear-weak',
                               nonzero.genes=nonzero.genes, n_batch=10, h0=26,
                               he_train=0,he_test=0,beta_sort_train=0,beta_sort_test=0,norm_type=2,
                               save_gene_data=T, save_surv_data=T, save_beta=T, plot_km=F, run_analysis=T,
                               seed=1234)
linear.weak.out[[3]]

# Squared terms
nl.qua.out = sim.survdata(surv.data,
                          N=10000, N.test=10000, test.ids=test.ids, sim.dataType='nl-quadratic', 
                          nonzero.genes=nonzero.genes, n_batch=10, h0=18, 
                          he_train=0, he_test=0, beta_sort_train=0, beta_sort_test=0, norm_type=2,
                          save_gene_data=T, save_surv_data=T, save_beta=T, plot_km=F, run_analysis=T,
                          seed=1234)
nl.qua.out[[3]]

# Squared terms with intercept
nl.shiftqua.out = sim.survdata(surv.data,
                               N=10000,N.test=10000,test.ids=test.ids,sim.dataType='nl-shiftquad',
                               nonzero.genes=nonzero.genes,n_batch=10,h0=21,
                               he_train=0,he_test=0,beta_sort_train=0,beta_sort_test=0,norm_type=2,
                               save_gene_data=F,save_surv_data=F,save_beta=F,plot_km=F,run_analysis=T,
                               seed=1234)
nl.shiftqua.out[[3]]
# Gene-gene interaction terms
nl.interact.out = sim.survdata(surv.data,
                               N=10000, N.test=10000, test.ids=test.ids, sim.dataType='nl-interaction', 
                               nonzero.genes=nonzero.genes, n_batch=10, h0=32, 
                               he_train=0, he_test=0, beta_sort_train=0, beta_sort_test=0, norm_type=2,
                               save_gene_data=T, save_surv_data=T, save_beta=T, plot_km=F, run_analysis=T,
                               seed=1234)
nl.interact.out[[3]]


# Sine interaction
nl.sine.out = sim.survdata(surv.data,
                           N=10000, N.test=10000, test.ids=test.ids, sim.dataType='nl-sine', 
                           nonzero.genes=nonzero.genes, n_batch=10, h0=24, 
                           he_train=0,he_test=0,beta_sort_train=0,beta_sort_test=0,norm_type=2,
                           save_gene_data=T, save_surv_data=T, save_beta=T, plot_km=F, run_analysis=T,
                           seed=1234)
nl.sine.out[[3]]


## ########################### BE00Asso00_normTMM ###############################
# Linear risk with moderate effects
linear.moderate.out = sim.survdata(surv.data,
                                   N=10000, N.test=10000, test.ids=test.ids, sim.dataType='linear-moderate',
                                   nonzero.genes=nonzero.genes, n_batch=10, h0=9.5,
                                   he_train=0,he_test=0,beta_sort_train=0,beta_sort_test=0, norm_type=3,
                                   save_gene_data=T, save_surv_data=T, save_beta=T, plot_km=F, run_analysis=T,
                                   1234)
linear.moderate.out[[3]] ## Look at the distribution of simulated data

#### NOTE: set save_gene_data = F from here
# Linear risk with weak effects
linear.weak.out = sim.survdata(surv.data,
                               N=10000, N.test=10000, test.ids=test.ids, sim.dataType='linear-weak',
                               nonzero.genes=nonzero.genes, n_batch=10, h0=26,
                               he_train=0,he_test=0,beta_sort_train=0,beta_sort_test=0,norm_type=3,
                               save_gene_data=T, save_surv_data=T, save_beta=T, plot_km=F, run_analysis=T,
                               seed=1234)
linear.weak.out[[3]]

# Squared terms
nl.qua.out = sim.survdata(surv.data,
                          N=10000, N.test=10000, test.ids=test.ids, sim.dataType='nl-quadratic', 
                          nonzero.genes=nonzero.genes, n_batch=10, h0=18, 
                          he_train=0, he_test=0, beta_sort_train=0, beta_sort_test=0, norm_type=3,
                          save_gene_data=T, save_surv_data=T, save_beta=T, plot_km=F, run_analysis=T,
                          seed=1234)
nl.qua.out[[3]]

# Squared terms with intercept
nl.shiftqua.out = sim.survdata(surv.data,
                               N=10000,N.test=10000,test.ids=test.ids,sim.dataType='nl-shiftquad',
                               nonzero.genes=nonzero.genes,n_batch=10,h0=21,
                               he_train=0,he_test=0,beta_sort_train=0,beta_sort_test=0,norm_type=3,
                               save_gene_data=F,save_surv_data=F,save_beta=F,plot_km=F,run_analysis=T,
                               seed=1234)
nl.shiftqua.out[[3]]

# Gene-gene interaction terms
nl.interact.out = sim.survdata(surv.data,
                               N=10000, N.test=10000, test.ids=test.ids, sim.dataType='nl-interaction', 
                               nonzero.genes=nonzero.genes, n_batch=10, h0=32, 
                               he_train=0, he_test=0, beta_sort_train=0, beta_sort_test=0, norm_type=3,
                               save_gene_data=T, save_surv_data=T, save_beta=T, plot_km=F, run_analysis=T,
                               seed=1234)
nl.interact.out[[3]]


# Sine interaction
nl.sine.out = sim.survdata(surv.data,
                           N=10000, N.test=10000, test.ids=test.ids, sim.dataType='nl-sine', 
                           nonzero.genes=nonzero.genes, n_batch=10, h0=24, 
                           he_train=0,he_test=0,beta_sort_train=0,beta_sort_test=0,norm_type=3,
                           save_gene_data=T, save_surv_data=T, save_beta=T, plot_km=F, run_analysis=T,
                           seed=1234)
nl.sine.out[[3]]


## ########################### BE00Asso00_normQuantile ###############################
# Linear risk with moderate effects
linear.moderate.out = sim.survdata(surv.data,
                                   N=10000, N.test=10000, test.ids=test.ids, sim.dataType='linear-moderate',
                                   nonzero.genes=nonzero.genes, n_batch=10, h0=9.5,
                                   he_train=0,he_test=0,beta_sort_train=0,beta_sort_test=0, norm_type=4,
                                   save_gene_data=F, save_surv_data=F, save_beta=F, plot_km=F,run_analysis=F,
                                   seed=1234)
linear.moderate.out[[3]] ## Look at the distribution of simulated data

#### NOTE: set save_gene_data = F from here
# Linear risk with weak effects
linear.weak.out = sim.survdata(surv.data,
                               N=10000, N.test=10000, test.ids=test.ids, sim.dataType='linear-weak',
                               nonzero.genes=nonzero.genes, n_batch=10, h0=26,
                               he_train=0,he_test=0,beta_sort_train=0,beta_sort_test=0,norm_type=4,
                               save_gene_data=F, save_surv_data=F, save_beta=F, plot_km=F,run_analysis=F,
                               seed=1234)
linear.weak.out[[3]]

# Squared terms
nl.qua.out = sim.survdata(surv.data,
                          N=10000, N.test=10000, test.ids=test.ids, sim.dataType='nl-quadratic', 
                          nonzero.genes=nonzero.genes, n_batch=10, h0=18, 
                          he_train=0, he_test=0, beta_sort_train=0, beta_sort_test=0, norm_type=4,
                          save_gene_data=F, save_surv_data=T, save_beta=T, plot_km=T,run_analysis=T,
                          seed=1234)
nl.qua.out[[3]]

# Squared terms with intercept
nl.shiftqua.out = sim.survdata(surv.data,
                               N=10000,N.test=10000,test.ids=test.ids,sim.dataType='nl-shiftquad',
                               nonzero.genes=nonzero.genes,n_batch=10,h0=21,
                               he_train=0,he_test=0,beta_sort_train=0,beta_sort_test=0,norm_type=4,
                               save_gene_data=T,save_surv_data=T,save_beta=T,plot_km=T,run_analysis=F,
                               seed=1234)
nl.shiftqua.out[[3]]

# Gene-gene interaction terms
nl.interact.out = sim.survdata(surv.data,
                               N=10000, N.test=10000, test.ids=test.ids, sim.dataType='nl-interaction', 
                               nonzero.genes=nonzero.genes, n_batch=10, h0=32, 
                               he_train=0, he_test=0, beta_sort_train=0, beta_sort_test=0, norm_type=4,
                               save_gene_data=F, save_surv_data=T, save_beta=T, plot_km=T,run_analysis=T,
                               seed=1234)
nl.interact.out[[3]]


# Sine interaction
nl.sine.out = sim.survdata(surv.data,
                           N=10000, N.test=10000, test.ids=test.ids, sim.dataType='nl-sine', 
                           nonzero.genes=nonzero.genes, n_batch=10, h0=24, 
                           he_train=0,he_test=0,beta_sort_train=0,beta_sort_test=0,norm_type=4,
                           save_gene_data=F, save_surv_data=T, save_beta=T, plot_km=T,run_analysis=T,
                           seed=1234)
nl.sine.out[[3]]



## ########################### BE00Asso00_normDEseq ###############################
# Linear risk with moderate effects
linear.moderate.out = sim.survdata(surv.data,
                                   N=10000, N.test=10000, test.ids=test.ids, sim.dataType='linear-moderate',
                                   nonzero.genes=nonzero.genes, n_batch=10, h0=13,
                                   he_train=0,he_test=0,beta_sort_train=0,beta_sort_test=0, norm_type=5,
                                   save_gene_data=F, save_surv_data=F, save_beta=F, plot_km=F,run_analysis=T,
                                   seed=1234)
linear.moderate.out[[3]] ## Look at the distribution of simulated data

#### NOTE: set save_gene_data = F from here
# Linear risk with weak effects
linear.weak.out = sim.survdata(surv.data,
                               N=10000, N.test=10000, test.ids=test.ids, sim.dataType='linear-weak',
                               nonzero.genes=nonzero.genes, n_batch=10, h0=26,
                               he_train=0,he_test=0,beta_sort_train=0,beta_sort_test=0,norm_type=5,
                               save_gene_data=F, save_surv_data=F, save_beta=F, plot_km=F,run_analysis=T,
                               seed=1234)
linear.weak.out[[3]]

# Squared terms
nl.qua.out = sim.survdata(surv.data,
                          N=10000, N.test=10000, test.ids=test.ids, sim.dataType='nl-quadratic', 
                          nonzero.genes=nonzero.genes, n_batch=10, h0=18, 
                          he_train=0, he_test=0, beta_sort_train=0, beta_sort_test=0, norm_type=5,
                          save_gene_data=F, save_surv_data=F, save_beta=F, plot_km=F,run_analysis=T,
                          seed=1234)
nl.qua.out[[3]]

# Squared terms with intercept
nl.shiftqua.out = sim.survdata(surv.data,
                               N=10000,N.test=10000,test.ids=test.ids,sim.dataType='nl-shiftquad',
                               nonzero.genes=nonzero.genes,n_batch=10,h0=21,
                               he_train=0,he_test=0,beta_sort_train=0,beta_sort_test=0,norm_type=5,
                               save_gene_data=T,save_surv_data=T,save_beta=T,plot_km=T,run_analysis=F,
                               seed=1234)
nl.shiftqua.out[[3]]
# Gene-gene interaction terms
nl.interact.out = sim.survdata(surv.data,
                               N=10000, N.test=10000, test.ids=test.ids, sim.dataType='nl-interaction', 
                               nonzero.genes=nonzero.genes, n_batch=10, h0=32, 
                               he_train=0, he_test=0, beta_sort_train=0, beta_sort_test=0, norm_type=5,
                               save_gene_data=F, save_surv_data=F, save_beta=F, plot_km=F,run_analysis=T,
                               seed=1234)
nl.interact.out[[3]]


# Sine interaction
nl.sine.out = sim.survdata(surv.data,
                           N=10000, N.test=10000, test.ids=test.ids, sim.dataType='nl-sine', 
                           nonzero.genes=nonzero.genes, n_batch=10, h0=24, 
                           he_train=0,he_test=0,beta_sort_train=0,beta_sort_test=0,norm_type=5,
                           save_gene_data=F, save_surv_data=F, save_beta=F, plot_km=F,run_analysis=T,
                           seed=1234)
nl.sine.out[[3]]






# ############################### BE10Asso00 #####################################
## ########################### BE10Asso00_normNone ###############################
# Linear risk with moderate effects
linear.moderate.out = sim.survdata(surv.data,
                                   N=10000, N.test=10000, test.ids=test.ids, sim.dataType='linear-moderate',
                                   nonzero.genes=nonzero.genes, n_batch=10, h0=13,
                                   he_train=1,he_test=0,beta_sort_train=0,beta_sort_test=0, norm_type=0,
                                   save_gene_data=T, save_surv_data=T, save_beta=T, plot_km=T,run_analysis=T,
                                   seed=1234)
linear.moderate.out[[3]] ## Look at the distribution of simulated data

#### NOTE: set save_gene_data = F from here
# Linear risk with weak effects
linear.weak.out = sim.survdata(surv.data,
                               N=10000, N.test=10000, test.ids=test.ids, sim.dataType='linear-weak',
                               nonzero.genes=nonzero.genes, n_batch=10, h0=26,
                               he_train=1,he_test=0,beta_sort_train=0,beta_sort_test=0,norm_type=0,
                               save_gene_data=F, save_surv_data=T, save_beta=T, plot_km=T,run_analysis=T,
                               seed=1234)
linear.weak.out[[3]]

# Squared terms
nl.qua.out = sim.survdata(surv.data,
                          N=10000, N.test=10000, test.ids=test.ids, sim.dataType='nl-quadratic', 
                          nonzero.genes=nonzero.genes, n_batch=10, h0=18, 
                          he_train=1, he_test=0, beta_sort_train=0, beta_sort_test=0, norm_type=0,
                          save_gene_data=F, save_surv_data=T, save_beta=T, plot_km=T,run_analysis=T,
                          seed=1234)
nl.qua.out[[3]]

# Squared terms with intercept
nl.shiftqua.out = sim.survdata(surv.data,
                               N=10000,N.test=10000,test.ids=test.ids,sim.dataType='nl-shiftquad',
                               nonzero.genes=nonzero.genes,n_batch=10,h0=21,
                               he_train=1,he_test=0,beta_sort_train=0,beta_sort_test=0,norm_type=0,
                               save_gene_data=T,save_surv_data=T,save_beta=T,plot_km=T,run_analysis=F,
                               seed=1234)
nl.shiftqua.out[[3]]

# Gene-gene interaction terms
nl.interact.out = sim.survdata(surv.data,
                               N=10000, N.test=10000, test.ids=test.ids, sim.dataType='nl-interaction', 
                               nonzero.genes=nonzero.genes, n_batch=10, h0=32, 
                               he_train=1, he_test=0, beta_sort_train=0, beta_sort_test=0, norm_type=0,
                               save_gene_data=F, save_surv_data=T, save_beta=T, plot_km=T,run_analysis=T,
                               seed=1234)
nl.interact.out[[3]]


# Sine interaction
nl.sine.out = sim.survdata(surv.data,
                           N=10000, N.test=10000, test.ids=test.ids, sim.dataType='nl-sine', 
                           nonzero.genes=nonzero.genes, n_batch=10, h0=24, 
                           he_train=1,he_test=0,beta_sort_train=0,beta_sort_test=0,norm_type=0,
                           save_gene_data=F, save_surv_data=T, save_beta=T, plot_km=T,run_analysis=T,
                           seed=1234)
nl.sine.out[[3]]



## ########################### BE10Asso00_normTC ###############################
# Linear risk with moderate effects
linear.moderate.out = sim.survdata(surv.data,
                                   N=10000, N.test=10000, test.ids=test.ids, sim.dataType='linear-moderate',
                                   nonzero.genes=nonzero.genes, n_batch=10, h0=13,
                                   he_train=1,he_test=0,beta_sort_train=0,beta_sort_test=0, norm_type=1,
                                   save_gene_data=T, save_surv_data=T, save_beta=T, plot_km=T,run_analysis=T,
                                   seed=1234)
linear.moderate.out[[3]] ## Look at the distribution of simulated data

#### NOTE: set save_gene_data = F from here
# Linear risk with weak effects
linear.weak.out = sim.survdata(surv.data,
                               N=10000, N.test=10000, test.ids=test.ids, sim.dataType='linear-weak',
                               nonzero.genes=nonzero.genes, n_batch=10, h0=26,
                               he_train=1,he_test=0,beta_sort_train=0,beta_sort_test=0,norm_type=1,
                               save_gene_data=F, save_surv_data=T, save_beta=T, plot_km=T,run_analysis=T,
                               seed=1234)
linear.weak.out[[3]]

# Squared terms
nl.qua.out = sim.survdata(surv.data,
                          N=10000, N.test=10000, test.ids=test.ids, sim.dataType='nl-quadratic', 
                          nonzero.genes=nonzero.genes, n_batch=10, h0=18, 
                          he_train=1, he_test=0, beta_sort_train=0, beta_sort_test=0, norm_type=1,
                          save_gene_data=F, save_surv_data=T, save_beta=T, plot_km=T,run_analysis=T,
                          seed=1234)
nl.qua.out[[3]]

# Squared terms with intercept
nl.shiftqua.out = sim.survdata(surv.data,
                               N=10000,N.test=10000,test.ids=test.ids,sim.dataType='nl-shiftquad',
                               nonzero.genes=nonzero.genes,n_batch=10,h0=21,
                               he_train=1,he_test=0,beta_sort_train=0,beta_sort_test=0,norm_type=1,
                               save_gene_data=T,save_surv_data=T,save_beta=T,plot_km=T,run_analysis=F,
                               seed=1234)
nl.shiftqua.out[[3]]

# Gene-gene interaction terms
nl.interact.out = sim.survdata(surv.data,
                               N=10000, N.test=10000, test.ids=test.ids, sim.dataType='nl-interaction', 
                               nonzero.genes=nonzero.genes, n_batch=10, h0=32, 
                               he_train=1, he_test=0, beta_sort_train=0, beta_sort_test=0, norm_type=1,
                               save_gene_data=F, save_surv_data=T, save_beta=T, plot_km=T,run_analysis=T,
                               seed=1234)
nl.interact.out[[3]]


# Sine interaction
nl.sine.out = sim.survdata(surv.data,
                           N=10000, N.test=10000, test.ids=test.ids, sim.dataType='nl-sine', 
                           nonzero.genes=nonzero.genes, n_batch=10, h0=24, 
                           he_train=1,he_test=0,beta_sort_train=0,beta_sort_test=0,norm_type=1,
                           save_gene_data=F, save_surv_data=T, save_beta=T, plot_km=T,run_analysis=T,
                           seed=1234)
nl.sine.out[[3]]


## ########################### BE10Asso00_normUQ ###############################
# Linear risk with moderate effects
linear.moderate.out = sim.survdata(surv.data,
                                   N=10000, N.test=10000, test.ids=test.ids, sim.dataType='linear-moderate',
                                   nonzero.genes=nonzero.genes, n_batch=10, h0=13,
                                   he_train=1,he_test=0,beta_sort_train=0,beta_sort_test=0, norm_type=2,
                                   save_gene_data=T, save_surv_data=T, save_beta=T, plot_km=T, run_analysis=T,
                                   seed=1234)
linear.moderate.out[[3]] ## Look at the distribution of simulated data

# Linear risk with weak effects
linear.weak.out = sim.survdata(surv.data,
                               N=10000, N.test=10000, test.ids=test.ids, sim.dataType='linear-weak',
                               nonzero.genes=nonzero.genes, n_batch=10, h0=26,
                               he_train=1,he_test=0,beta_sort_train=0,beta_sort_test=0,norm_type=2,
                               save_gene_data=F, save_surv_data=T, save_beta=T, plot_km=T,run_analysis=T,
                               seed=1234)
linear.weak.out[[3]]

# Squared terms
nl.qua.out = sim.survdata(surv.data,
                          N=10000, N.test=10000, test.ids=test.ids, sim.dataType='nl-quadratic',
                          nonzero.genes=nonzero.genes, n_batch=10, h0=18,
                          he_train=1, he_test=0, beta_sort_train=0, beta_sort_test=0, norm_type=2,
                          save_gene_data=F, save_surv_data=T, save_beta=T, plot_km=T,run_analysis=T,
                          seed=1234)
nl.qua.out[[3]]

# Squared terms with intercept
nl.shiftqua.out = sim.survdata(surv.data,
                               N=10000,N.test=10000,test.ids=test.ids,sim.dataType='nl-shiftquad',
                               nonzero.genes=nonzero.genes,n_batch=10,h0=21,
                               he_train=1,he_test=0,beta_sort_train=0,beta_sort_test=0,norm_type=2,
                               save_gene_data=T,save_surv_data=T,save_beta=T,plot_km=T,run_analysis=F,
                               seed=1234)
nl.shiftqua.out[[3]]

# Gene-gene interaction terms
nl.interact.out = sim.survdata(surv.data,
                               N=10000, N.test=10000, test.ids=test.ids, sim.dataType='nl-interaction', 
                               nonzero.genes=nonzero.genes, n_batch=10, h0=32,
                               he_train=1, he_test=0, beta_sort_train=0, beta_sort_test=0, norm_type=2,
                               save_gene_data=F, save_surv_data=T, save_beta=T, plot_km=T,run_analysis=T,
                               seed=1234)
nl.interact.out[[3]]

# Sine interaction
nl.sine.out = sim.survdata(surv.data,
                           N=10000, N.test=10000, test.ids=test.ids, sim.dataType='nl-sine',
                           nonzero.genes=nonzero.genes, n_batch=10, h0=24,
                           he_train=1,he_test=0,beta_sort_train=0,beta_sort_test=0,norm_type=2,
                           save_gene_data=F, save_surv_data=T, save_beta=T, plot_km=T,run_analysis=T,
                           seed=1234)
nl.sine.out[[3]]



## ########################### BE10Asso00_normTMM ###############################
# Linear risk with moderate effects
# linear.moderate.out = sim.survdata(surv.data,
#                                    N=10000, N.test=10000, test.ids=test.ids, sim.dataType='linear-moderate',
#                                    nonzero.genes=nonzero.genes, n_batch=10, h0=13,
#                                    he_train=1,he_test=0,beta_sort_train=0,beta_sort_test=0, norm_type=3,
#                                    save_gene_data=T, save_surv_data=T, save_beta=T, plot_km=T, run_analysis=T,
#                                    seed=1234)
# linear.moderate.out[[3]] ## Look at the distribution of simulated data

# # Linear risk with weak effects
# linear.weak.out = sim.survdata(surv.data,
#                                N=10000, N.test=10000, test.ids=test.ids, sim.dataType='linear-weak',
#                                nonzero.genes=nonzero.genes, n_batch=10, h0=26,
#                                he_train=1,he_test=0,beta_sort_train=0,beta_sort_test=0,norm_type=3,
#                                save_gene_data=F, save_surv_data=T, save_beta=T, plot_km=T,run_analysis=T,
#                                seed=1234)
# linear.weak.out[[3]]
# 
# # Squared terms
# nl.qua.out = sim.survdata(surv.data,
#                           N=10000, N.test=10000, test.ids=test.ids, sim.dataType='nl-quadratic', 
#                           nonzero.genes=nonzero.genes, n_batch=10, h0=18, 
#                           he_train=1, he_test=0, beta_sort_train=0, beta_sort_test=0, norm_type=3,
#                           save_gene_data=F, save_surv_data=T, save_beta=T, plot_km=T,run_analysis=T,
#                           seed=1234)
# nl.qua.out[[3]]

# Squared terms with intercept
nl.shiftqua.out = sim.survdata(surv.data,
                               N=10000,N.test=10000,test.ids=test.ids,sim.dataType='nl-shiftquad',
                               nonzero.genes=nonzero.genes,n_batch=10,h0=21,
                               he_train=1,he_test=0,beta_sort_train=0,beta_sort_test=0,norm_type=3,
                               save_gene_data=T,save_surv_data=T,save_beta=T,plot_km=T,run_analysis=F,
                               seed=1234)
nl.shiftqua.out[[3]]

# Gene-gene interaction terms
# nl.interact.out = sim.survdata(surv.data,
#                                N=10000, N.test=10000, test.ids=test.ids, sim.dataType='nl-interaction', 
#                                nonzero.genes=nonzero.genes, n_batch=10, h0=32, 
#                                he_train=1, he_test=0, beta_sort_train=0, beta_sort_test=0, norm_type=3,
#                                save_gene_data=F, save_surv_data=T, save_beta=T, plot_km=T,run_analysis=T,
#                                seed=1234)
# nl.interact.out[[3]]
# 
# # Sine interaction
# nl.sine.out = sim.survdata(surv.data,
#                            N=10000, N.test=10000, test.ids=test.ids, sim.dataType='nl-sine', 
#                            nonzero.genes=nonzero.genes, n_batch=10, h0=24, 
#                            he_train=1,he_test=0,beta_sort_train=0,beta_sort_test=0,norm_type=3,
#                            save_gene_data=F, save_surv_data=T, save_beta=T, plot_km=T,run_analysis=T,
#                            seed=1234)
# nl.sine.out[[3]]



## ########################### BE10Asso00_normQuantile ###############################
# Linear risk with moderate effects
linear.moderate.out = sim.survdata(surv.data,
                                   N=10000, N.test=10000, test.ids=test.ids, sim.dataType='linear-moderate',
                                   nonzero.genes=nonzero.genes, n_batch=10, h0=13,
                                   he_train=1,he_test=0,beta_sort_train=0,beta_sort_test=0, norm_type=4,
                                   save_gene_data=T, save_surv_data=T, save_beta=T, plot_km=T, run_analysis=T,
                                   seed=1234)
linear.moderate.out[[3]] ## Look at the distribution of simulated data

# Linear risk with weak effects
linear.weak.out = sim.survdata(surv.data,
                               N=10000, N.test=10000, test.ids=test.ids, sim.dataType='linear-weak',
                               nonzero.genes=nonzero.genes, n_batch=10, h0=26,
                               he_train=1,he_test=0,beta_sort_train=0,beta_sort_test=0,norm_type=4,
                               save_gene_data=F, save_surv_data=T, save_beta=T, plot_km=T,run_analysis=T,
                               seed=1234)
linear.weak.out[[3]]

# Squared terms
nl.qua.out = sim.survdata(surv.data,
                          N=10000, N.test=10000, test.ids=test.ids, sim.dataType='nl-quadratic', 
                          nonzero.genes=nonzero.genes, n_batch=10, h0=18, 
                          he_train=1, he_test=0, beta_sort_train=0, beta_sort_test=0, norm_type=4,
                          save_gene_data=F, save_surv_data=T, save_beta=T, plot_km=T,run_analysis=T,
                          seed=1234)
nl.qua.out[[3]]

# Squared terms with intercept
nl.shiftqua.out = sim.survdata(surv.data,
                               N=10000,N.test=10000,test.ids=test.ids,sim.dataType='nl-shiftquad',
                               nonzero.genes=nonzero.genes,n_batch=10,h0=21,
                               he_train=1,he_test=0,beta_sort_train=0,beta_sort_test=0,norm_type=4,
                               save_gene_data=T,save_surv_data=T,save_beta=T,plot_km=T,run_analysis=F,
                               seed=1234)
nl.shiftqua.out[[3]]

# Gene-gene interaction terms
nl.interact.out = sim.survdata(surv.data,
                               N=10000, N.test=10000, test.ids=test.ids, sim.dataType='nl-interaction', 
                               nonzero.genes=nonzero.genes, n_batch=10, h0=32, 
                               he_train=1, he_test=0, beta_sort_train=0, beta_sort_test=0, norm_type=4,
                               save_gene_data=F, save_surv_data=T, save_beta=T, plot_km=T,run_analysis=T,
                               seed=1234)
nl.interact.out[[3]]

# Sine interaction
nl.sine.out = sim.survdata(surv.data,
                           N=10000, N.test=10000, test.ids=test.ids, sim.dataType='nl-sine', 
                           nonzero.genes=nonzero.genes, n_batch=10, h0=24, 
                           he_train=1,he_test=0,beta_sort_train=0,beta_sort_test=0,norm_type=4,
                           save_gene_data=F, save_surv_data=T, save_beta=T, plot_km=T,run_analysis=T,
                           seed=1234)
nl.sine.out[[3]]



## ########################### BE10Asso00_normDEseq ###############################
# Linear risk with moderate effects
linear.moderate.out = sim.survdata(surv.data,
                                   N=10000, N.test=10000, test.ids=test.ids, sim.dataType='linear-moderate',
                                   nonzero.genes=nonzero.genes, n_batch=10, h0=13,
                                   he_train=1,he_test=0,beta_sort_train=0,beta_sort_test=0, norm_type=5,
                                   save_gene_data=T, save_surv_data=T, save_beta=T, plot_km=T, run_analysis=T,
                                   seed=1234)
linear.moderate.out[[3]] ## Look at the distribution of simulated data

# Linear risk with weak effects
linear.weak.out = sim.survdata(surv.data,
                               N=10000, N.test=10000, test.ids=test.ids, sim.dataType='linear-weak',
                               nonzero.genes=nonzero.genes, n_batch=10, h0=26,
                               he_train=1,he_test=0,beta_sort_train=0,beta_sort_test=0,norm_type=5,
                               save_gene_data=F, save_surv_data=T, save_beta=T, plot_km=T,run_analysis=T,
                               seed=1234)
linear.weak.out[[3]]

# Squared terms
nl.qua.out = sim.survdata(surv.data,
                          N=10000, N.test=10000, test.ids=test.ids, sim.dataType='nl-quadratic', 
                          nonzero.genes=nonzero.genes, n_batch=10, h0=18, 
                          he_train=1, he_test=0, beta_sort_train=0, beta_sort_test=0, norm_type=5,
                          save_gene_data=F, save_surv_data=T, save_beta=T, plot_km=T,run_analysis=T,
                          seed=1234)
nl.qua.out[[3]]

# Squared terms with intercept
nl.shiftqua.out = sim.survdata(surv.data,
                               N=10000,N.test=10000,test.ids=test.ids,sim.dataType='nl-shiftquad',
                               nonzero.genes=nonzero.genes,n_batch=10,h0=21,
                               he_train=1,he_test=0,beta_sort_train=0,beta_sort_test=0,norm_type=5,
                               save_gene_data=T,save_surv_data=T,save_beta=T,plot_km=T,run_analysis=F,
                               seed=1234)
nl.shiftqua.out[[3]]

# Gene-gene interaction terms
nl.interact.out = sim.survdata(surv.data,
                               N=10000, N.test=10000, test.ids=test.ids, sim.dataType='nl-interaction', 
                               nonzero.genes=nonzero.genes, n_batch=10, h0=32, 
                               he_train=1, he_test=0, beta_sort_train=0, beta_sort_test=0, norm_type=5,
                               save_gene_data=F, save_surv_data=T, save_beta=T, plot_km=T,run_analysis=T,
                               seed=1234)
nl.interact.out[[3]]

# Sine interaction
nl.sine.out = sim.survdata(surv.data,
                           N=10000, N.test=10000, test.ids=test.ids, sim.dataType='nl-sine', 
                           nonzero.genes=nonzero.genes, n_batch=10, h0=24, 
                           he_train=1,he_test=0,beta_sort_train=0,beta_sort_test=0,norm_type=5,
                           save_gene_data=F, save_surv_data=T, save_beta=T, plot_km=T,run_analysis=T,
                           seed=1234)
nl.sine.out[[3]]





# ############################### BE11Asso00 #####################################
## ########################### BE11Asso00_normNone ###############################
# Linear risk with moderate effects
linear.moderate.out = sim.survdata(surv.data,
                                   N=10000,N.test=10000,test.ids=test.ids,sim.dataType='linear-moderate',
                                   nonzero.genes=nonzero.genes,n_batch=10,h0=13,
                                   he_train=1,he_test=1,beta_sort_train=0,beta_sort_test=0,norm_type=0,
                                   save_gene_data=T,save_surv_data=T,save_beta=T,plot_km=T,run_analysis=T,
                                   seed=1234)
linear.moderate.out[[3]] ## Look at the distribution of simulated data

#### NOTE: set save_gene_data = F from here
# Linear risk with weak effects
linear.weak.out = sim.survdata(surv.data,
                               N=10000,N.test=10000,test.ids=test.ids,sim.dataType='linear-weak',
                               nonzero.genes=nonzero.genes,n_batch=10,h0=26,
                               he_train=1,he_test=1,beta_sort_train=0,beta_sort_test=0,norm_type=0,
                               save_gene_data=F,save_surv_data=T,save_beta=T,plot_km=T,run_analysis=T,
                               seed=1234)
linear.weak.out[[3]]

# Squared terms
nl.qua.out = sim.survdata(surv.data,
                          N=10000,N.test=10000,test.ids=test.ids,sim.dataType='nl-quadratic',
                          nonzero.genes=nonzero.genes,n_batch=10,h0=18,
                          he_train=1,he_test=1,beta_sort_train=0,beta_sort_test=0,norm_type=0,
                          save_gene_data=F,save_surv_data=T,save_beta=T,plot_km=T,run_analysis=T,
                          seed=1234)
nl.qua.out[[3]]


# Gene-gene interaction terms
nl.interact.out = sim.survdata(surv.data,
                               N=10000,N.test=10000,test.ids=test.ids,sim.dataType='nl-interaction',
                               nonzero.genes=nonzero.genes,n_batch=10,h0=32,
                               he_train=1,he_test=1,beta_sort_train=0,beta_sort_test=0,norm_type=0,
                               save_gene_data=F,save_surv_data=T,save_beta=T,plot_km=T,run_analysis=T,
                               seed=1234)
nl.interact.out[[3]]


# Sine interaction
nl.sine.out = sim.survdata(surv.data,
                           N=10000,N.test=10000,test.ids=test.ids,sim.dataType='nl-sine',
                           nonzero.genes=nonzero.genes,n_batch=10,h0=24,
                           he_train=1,he_test=1,beta_sort_train=0,beta_sort_test=0,norm_type=0,
                           save_gene_data=F,save_surv_data=T,save_beta=T,plot_km=T,run_analysis=T,
                           seed=1234)
nl.sine.out[[3]]



## ########################### BE11Asso00_normTC ###############################
# Linear risk with moderate effects
linear.moderate.out = sim.survdata(surv.data,
                                   N=10000,N.test=10000,test.ids=test.ids,sim.dataType='linear-moderate',
                                   nonzero.genes=nonzero.genes,n_batch=10,h0=13,
                                   he_train=1,he_test=1,beta_sort_train=0,beta_sort_test=0,norm_type=1,
                                   save_gene_data=T,save_surv_data=T,save_beta=T,plot_km=T,run_analysis=T,
                                   seed=1234)
linear.moderate.out[[3]] ## Look at the distribution of simulated data

#### NOTE: set save_gene_data = F from here
# Linear risk with weak effects
linear.weak.out = sim.survdata(surv.data,
                               N=10000,N.test=10000,test.ids=test.ids,sim.dataType='linear-weak',
                               nonzero.genes=nonzero.genes,n_batch=10,h0=26,
                               he_train=1,he_test=1,beta_sort_train=0,beta_sort_test=0,norm_type=1,
                               save_gene_data=F,save_surv_data=T,save_beta=T,plot_km=T,run_analysis=T,
                               seed=1234)
linear.weak.out[[3]]

# Squared terms
nl.qua.out = sim.survdata(surv.data,
                          N=10000,N.test=10000,test.ids=test.ids,sim.dataType='nl-quadratic',
                          nonzero.genes=nonzero.genes,n_batch=10,h0=18,
                          he_train=1,he_test=1,beta_sort_train=0,beta_sort_test=0,norm_type=1,
                          save_gene_data=F,save_surv_data=T,save_beta=T,plot_km=T,run_analysis=T,
                          seed=1234)
nl.qua.out[[3]]


# Gene-gene interaction terms
nl.interact.out = sim.survdata(surv.data,
                               N=10000,N.test=10000,test.ids=test.ids,sim.dataType='nl-interaction',
                               nonzero.genes=nonzero.genes,n_batch=10,h0=32,
                               he_train=1,he_test=1,beta_sort_train=0,beta_sort_test=0,norm_type=1,
                               save_gene_data=F,save_surv_data=T,save_beta=T,plot_km=T,run_analysis=T,
                               seed=1234)
nl.interact.out[[3]]


# Sine interaction
nl.sine.out = sim.survdata(surv.data,
                           N=10000,N.test=10000,test.ids=test.ids,sim.dataType='nl-sine',
                           nonzero.genes=nonzero.genes,n_batch=10,h0=24,
                           he_train=1,he_test=1,beta_sort_train=0,beta_sort_test=0,norm_type=1,
                           save_gene_data=F,save_surv_data=T,save_beta=T,plot_km=T,run_analysis=T,
                           seed=1234)
nl.sine.out[[3]]


## ########################### BE11Asso00_normUQ ###############################
# Linear risk with moderate effects
linear.moderate.out = sim.survdata(surv.data,
                                   N=10000,N.test=10000,test.ids=test.ids,sim.dataType='linear-moderate',
                                   nonzero.genes=nonzero.genes,n_batch=10,h0=13,
                                   he_train=1,he_test=1,beta_sort_train=0,beta_sort_test=0,norm_type=2,
                                   save_gene_data=T,save_surv_data=T,save_beta=T,plot_km=T,run_analysis=F,
                                   seed=1234)
linear.moderate.out[[3]] ## Look at the distribution of simulated data

# Linear risk with weak effects
linear.weak.out = sim.survdata(surv.data,
                               N=10000,N.test=10000,test.ids=test.ids,sim.dataType='linear-weak',
                               nonzero.genes=nonzero.genes,n_batch=10,h0=26,
                               he_train=1,he_test=1,beta_sort_train=0,beta_sort_test=0,norm_type=2,
                               save_gene_data=F,save_surv_data=T,save_beta=T,plot_km=T,run_analysis=F,
                               seed=1234)
linear.weak.out[[3]]

# Squared terms
nl.qua.out = sim.survdata(surv.data,
                          N=10000,N.test=10000,test.ids=test.ids,sim.dataType='nl-quadratic',
                          nonzero.genes=nonzero.genes,n_batch=10,h0=18,
                          he_train=1,he_test=1,beta_sort_train=0,beta_sort_test=0,norm_type=2,
                          save_gene_data=F,save_surv_data=T,save_beta=T,plot_km=T,run_analysis=F,
                          seed=1234)
nl.qua.out[[3]]

# Gene-gene interaction terms
nl.interact.out = sim.survdata(surv.data,
                               N=10000,N.test=10000,test.ids=test.ids,sim.dataType='nl-interaction',
                               nonzero.genes=nonzero.genes,n_batch=10,h0=32,
                               he_train=1,he_test=1,beta_sort_train=0,beta_sort_test=0,norm_type=2,
                               save_gene_data=F,save_surv_data=T,save_beta=T,plot_km=T,run_analysis=F,
                               seed=1234)
nl.interact.out[[3]]

# Sine interaction
nl.sine.out = sim.survdata(surv.data,
                           N=10000,N.test=10000,test.ids=test.ids,sim.dataType='nl-sine',
                           nonzero.genes=nonzero.genes,n_batch=10,h0=24,
                           he_train=1,he_test=1,beta_sort_train=0,beta_sort_test=0,norm_type=2,
                           save_gene_data=F,save_surv_data=T,save_beta=T,plot_km=T,run_analysis=F,
                           seed=1234)
nl.sine.out[[3]]



## ########################### BE11Asso00_normTMM ###############################
# Linear risk with moderate effects
linear.moderate.out = sim.survdata(surv.data,
                                   N=10000,N.test=10000,test.ids=test.ids,sim.dataType='linear-moderate',
                                   nonzero.genes=nonzero.genes,n_batch=10,h0=13,
                                   he_train=1,he_test=1,beta_sort_train=0,beta_sort_test=0,norm_type=3,
                                   save_gene_data=T,save_surv_data=T,save_beta=T,plot_km=T,run_analysis=F,
                                   seed=1234)
linear.moderate.out[[3]] ## Look at the distribution of simulated data

# Linear risk with weak effects
linear.weak.out = sim.survdata(surv.data,
                               N=10000,N.test=10000,test.ids=test.ids,sim.dataType='linear-weak',
                               nonzero.genes=nonzero.genes,n_batch=10,h0=26,
                               he_train=1,he_test=1,beta_sort_train=0,beta_sort_test=0,norm_type=3,
                               save_gene_data=F,save_surv_data=T,save_beta=T,plot_km=T,run_analysis=F,
                               seed=1234)
linear.weak.out[[3]]

# Squared terms
nl.qua.out = sim.survdata(surv.data,
                          N=10000,N.test=10000,test.ids=test.ids,sim.dataType='nl-quadratic',
                          nonzero.genes=nonzero.genes,n_batch=10,h0=18,
                          he_train=1,he_test=1,beta_sort_train=0,beta_sort_test=0,norm_type=3,
                          save_gene_data=F,save_surv_data=T,save_beta=T,plot_km=T,run_analysis=F,
                          seed=1234)
nl.qua.out[[3]]

# Squared terms with intercept
nl.qua.out = sim.survdata(surv.data,
                          N=10000,N.test=10000,test.ids=test.ids,sim.dataType='nl-shiftquad',
                          nonzero.genes=nonzero.genes,n_batch=10,h0=18,
                          he_train=1,he_test=1,beta_sort_train=0,beta_sort_test=0,norm_type=3,
                          save_gene_data=F,save_surv_data=T,save_beta=T,plot_km=T,run_analysis=F,
                          seed=1234)
nl.qua.out[[3]]

# Gene-gene interaction terms
nl.interact.out = sim.survdata(surv.data,
                               N=10000,N.test=10000,test.ids=test.ids,sim.dataType='nl-interaction',
                               nonzero.genes=nonzero.genes,n_batch=10,h0=32,
                               he_train=1,he_test=1,beta_sort_train=0,beta_sort_test=0,norm_type=3,
                               save_gene_data=F,save_surv_data=T,save_beta=T,plot_km=T,run_analysis=F,
                               seed=1234)
nl.interact.out[[3]]

# Sine interaction
nl.sine.out = sim.survdata(surv.data,
                           N=10000,N.test=10000,test.ids=test.ids,sim.dataType='nl-sine',
                           nonzero.genes=nonzero.genes,n_batch=10,h0=24,
                           he_train=1,he_test=1,beta_sort_train=0,beta_sort_test=0,norm_type=3,
                           save_gene_data=F,save_surv_data=T,save_beta=T,plot_km=T,run_analysis=F,
                           seed=1234)
nl.sine.out[[3]]



## ########################### BE11Asso00_normQuantile ###############################
# Linear risk with moderate effects
linear.moderate.out = sim.survdata(surv.data,
                                   N=10000,N.test=10000,test.ids=test.ids,sim.dataType='linear-moderate',
                                   nonzero.genes=nonzero.genes,n_batch=10,h0=13,
                                   he_train=1,he_test=1,beta_sort_train=0,beta_sort_test=0,norm_type=4,
                                   save_gene_data=T,save_surv_data=T,save_beta=T,plot_km=T,run_analysis=F,
                                   seed=1234)
linear.moderate.out[[3]] ## Look at the distribution of simulated data

# Linear risk with weak effects
linear.weak.out = sim.survdata(surv.data,
                               N=10000,N.test=10000,test.ids=test.ids,sim.dataType='linear-weak',
                               nonzero.genes=nonzero.genes,n_batch=10,h0=26,
                               he_train=1,he_test=1,beta_sort_train=0,beta_sort_test=0,norm_type=4,
                               save_gene_data=F,save_surv_data=T,save_beta=T,plot_km=T,run_analysis=F,
                               seed=1234)
linear.weak.out[[3]]

# Squared terms
nl.qua.out = sim.survdata(surv.data,
                          N=10000,N.test=10000,test.ids=test.ids,sim.dataType='nl-quadratic',
                          nonzero.genes=nonzero.genes,n_batch=10,h0=18,
                          he_train=1,he_test=1,beta_sort_train=0,beta_sort_test=0,norm_type=4,
                          save_gene_data=F,save_surv_data=T,save_beta=T,plot_km=T,run_analysis=F,
                          seed=1234)
nl.qua.out[[3]]

# Gene-gene interaction terms
nl.interact.out = sim.survdata(surv.data,
                               N=10000,N.test=10000,test.ids=test.ids,sim.dataType='nl-interaction',
                               nonzero.genes=nonzero.genes,n_batch=10,h0=32,
                               he_train=1,he_test=1,beta_sort_train=0,beta_sort_test=0,norm_type=4,
                               save_gene_data=F,save_surv_data=T,save_beta=T,plot_km=T,run_analysis=F,
                               seed=1234)
nl.interact.out[[3]]

# Sine interaction
nl.sine.out = sim.survdata(surv.data,
                           N=10000,N.test=10000,test.ids=test.ids,sim.dataType='nl-sine',
                           nonzero.genes=nonzero.genes,n_batch=10,h0=24,
                           he_train=1,he_test=1,beta_sort_train=0,beta_sort_test=0,norm_type=4,
                           save_gene_data=F,save_surv_data=T,save_beta=T,plot_km=T,run_analysis=F,
                           seed=1234)
nl.sine.out[[3]]



## ########################### BE11Asso00_normDEseq ###############################
# Linear risk with moderate effects
linear.moderate.out = sim.survdata(surv.data,
                                   N=10000,N.test=10000,test.ids=test.ids,sim.dataType='linear-moderate',
                                   nonzero.genes=nonzero.genes,n_batch=10,h0=13,
                                   he_train=1,he_test=1,beta_sort_train=0,beta_sort_test=0,norm_type=5,
                                   save_gene_data=T,save_surv_data=T,save_beta=T,plot_km=T,run_analysis=F,
                                   seed=1234)
linear.moderate.out[[3]] ## Look at the distribution of simulated data

# Linear risk with weak effects
linear.weak.out = sim.survdata(surv.data,
                               N=10000,N.test=10000,test.ids=test.ids,sim.dataType='linear-weak',
                               nonzero.genes=nonzero.genes,n_batch=10,h0=26,
                               he_train=1,he_test=1,beta_sort_train=0,beta_sort_test=0,norm_type=5,
                               save_gene_data=F,save_surv_data=T,save_beta=T,plot_km=T,run_analysis=F,
                               seed=1234)
linear.weak.out[[3]]

# Squared terms
nl.qua.out = sim.survdata(surv.data,
                          N=10000,N.test=10000,test.ids=test.ids,sim.dataType='nl-quadratic',
                          nonzero.genes=nonzero.genes,n_batch=10,h0=18,
                          he_train=1,he_test=1,beta_sort_train=0,beta_sort_test=0,norm_type=5,
                          save_gene_data=F,save_surv_data=T,save_beta=T,plot_km=T,run_analysis=F,
                          seed=1234)
nl.qua.out[[3]]

# Gene-gene interaction terms
nl.interact.out = sim.survdata(surv.data,
                               N=10000,N.test=10000,test.ids=test.ids,sim.dataType='nl-interaction',
                               nonzero.genes=nonzero.genes,n_batch=10,h0=32,
                               he_train=1,he_test=1,beta_sort_train=0,beta_sort_test=0,norm_type=5,
                               save_gene_data=F,save_surv_data=T,save_beta=T,plot_km=T,run_analysis=F,
                               seed=1234)
nl.interact.out[[3]]

# Sine interaction
nl.sine.out = sim.survdata(surv.data,
                           N=10000,N.test=10000,test.ids=test.ids,sim.dataType='nl-sine',
                           nonzero.genes=nonzero.genes,n_batch=10,h0=24,
                           he_train=1,he_test=1,beta_sort_train=0,beta_sort_test=0,norm_type=5,
                           save_gene_data=F,save_surv_data=T,save_beta=T,plot_km=T,run_analysis=F,
                           seed=1234)
nl.sine.out[[3]]
