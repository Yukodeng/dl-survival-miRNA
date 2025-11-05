# ########################## Run Simulation ##################################
# options(echo=TRUE)
# NOTE: supply job array number 1-7
# norm_type <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID')) - 1

source('~/dl-survival-miRNA/simulation/aug-batch-survival-pipeline-v2.R')


## ############################# Test ###############################
# h0.mod <- 2.75
# h0.weak <- 65
# h0.qua <- 22
# h0.shift <- 19
# h0.inter <- 25
# h0.sine <- 0.8
# norm_type = 0
# # Linear risk with moderate effects - BE10Asso10_normNone
# linear.moderate.out = sim.survdata(surv.data,
#                                    N=10000, N.test=10000, test.ids=test.ids, sim.dataType='linear-moderate',
#                                    nonzero.genes=nonzero.genes, n_batch=10, h0=2.75,
#                                    he_train=1,he_test=0,beta_sort_train=0.05,beta_sort_test=0,norm_type=norm_type,
#                                    save_gene_data=F, save_surv_data=F, save_beta=F, plot_km=F, run_analysis=T,
#                                    subset_sizes = c(100,5000),
#                                    runs_per_size = c(20,10),
#                                    splits_per_size = c(3,10),
#                                    seed=1234)


# ## ############################# BE00Asso00 ###############################
# for (norm_type in 0:6) {
# 
#   h0.mod <- 2.75
#   h0.weak <- 65
#   h0.qua <- 22
#   h0.shift <- 19
#   h0.inter <- 25
#   h0.sine <- 0.8
# 
#   if (norm_type==5) {
#     h0.weak = 9.5
#     h0.qua = 15
#     h0.shift = 16
#   }
# 
#   
#     # Linear risk with moderate effects
#     linear.moderate.out = sim.survdata(surv.data,
#                                        N=10000, N.test=10000, test.ids=test.ids, sim.dataType='linear-moderate',
#                                        nonzero.genes=nonzero.genes, n_batch=10, h0=h0.mod,
#                                        he_train=0,he_test=0,beta_sort_train=0,beta_sort_test=0,norm_type=norm_type,
#                                        save_gene_data=F, save_surv_data=F, save_beta=F, plot_km=F, run_analysis=T,
#                                        seed=1234)
#     # Linear risk with weak effects
#     linear.weak.out = sim.survdata(surv.data,
#                                    N=10000, N.test=10000, test.ids=test.ids, sim.dataType='linear-weak',
#                                    nonzero.genes=nonzero.genes, n_batch=10, h0=h0.weak,
#                                    he_train=0,he_test=0,beta_sort_train=0,beta_sort_test=0,norm_type=norm_type,
#                                    save_gene_data=F,save_surv_data=F, save_beta=F, plot_km=F, run_analysis=T,
#                                    seed=1234)
#     # Squared terms
#     nl.qua.out = sim.survdata(surv.data,
#                               N=10000, N.test=10000, test.ids=test.ids, sim.dataType='nl-quadratic',
#                               nonzero.genes=nonzero.genes, n_batch=10, h0=h0.qua,
#                               he_train=0,he_test=0,beta_sort_train=0,beta_sort_test=0,norm_type=norm_type,
#                               save_gene_data=F, save_surv_data=F, save_beta=F, plot_km=F, run_analysis=T,
#                               seed=1234)
#     
#     # Squared terms with intercept
#     nl.shiftqua.out = sim.survdata(surv.data,
#                                    N=10000,N.test=10000,test.ids=test.ids,sim.dataType='nl-shiftquad',
#                                    nonzero.genes=nonzero.genes,n_batch=10,h0=h0.shift,
#                                    he_train=0,he_test=0,beta_sort_train=0,beta_sort_test=0,norm_type=norm_type,
#                                    save_gene_data=F, save_surv_data=F, save_beta=F, plot_km=F, run_analysis=T,
#                                    seed=1234)
#   
#   # Gene-gene interaction terms
#   nl.interact.out = sim.survdata(surv.data,
#                                  N=10000, N.test=10000, test.ids=test.ids, sim.dataType='nl-interaction',
#                                  nonzero.genes=nonzero.genes, n_batch=10, h0=h0.inter,
#                                  he_train=0,he_test=0,beta_sort_train=0,beta_sort_test=0,norm_type=norm_type,
#                                  save_gene_data=F, save_surv_data=F, save_beta=F, plot_km=F, run_analysis=T,
#                                  seed=1234)
#   # Sine interaction
#   nl.sine.out = sim.survdata(surv.data,
#                              N=10000, N.test=10000, test.ids=test.ids, sim.dataType='nl-sine',
#                              nonzero.genes=nonzero.genes, n_batch=10, h0=h0.sine,
#                              he_train=0,he_test=0,beta_sort_train=0,beta_sort_test=0,norm_type=norm_type,
#                              save_gene_data=F, save_surv_data=F, save_beta=F, plot_km=F, run_analysis=T,
#                              seed=1234)
# }


# ## ############################# BE10Asso00 ####################################
# for (norm_type in 0:6) {
# 
#   h0.mod <- 2.75
#   h0.weak <- 65
#   h0.qua <- 22
#   h0.shift <- 19
#   h0.inter <- 25
#   h0.sine <- 0.8
# 
#   if (norm_type==5) {
#     h0.weak = 9.5
#     h0.qua = 15
#     h0.shift = 16
#   }
# 
#   # Linear risk with moderate effects
#   linear.moderate.out = sim.survdata(surv.data,
#                                      N=10000, N.test=10000, test.ids=test.ids, sim.dataType='linear-moderate',
#                                      nonzero.genes=nonzero.genes, n_batch=10, h0=h0.mod,
#                                      he_train=1,he_test=0,beta_sort_train=0,beta_sort_test=0,norm_type=norm_type,
#                                      save_gene_data=F, save_surv_data=F, save_beta=F, plot_km=F, run_analysis=T,
#                                      seed=1234)
#   # Linear risk with weak effects
#   linear.weak.out = sim.survdata(surv.data,
#                                  N=10000, N.test=10000, test.ids=test.ids, sim.dataType='linear-weak',
#                                  nonzero.genes=nonzero.genes, n_batch=10, h0=h0.weak,
#                                  he_train=1,he_test=0,beta_sort_train=0,beta_sort_test=0,norm_type=norm_type,
#                                  save_gene_data=F, save_surv_data=F, save_beta=F, plot_km=F, run_analysis=T,
#                                  seed=1234)
#   # Squared terms
#   nl.qua.out = sim.survdata(surv.data,
#                             N=10000, N.test=10000, test.ids=test.ids, sim.dataType='nl-quadratic',
#                             nonzero.genes=nonzero.genes, n_batch=10, h0=h0.qua,
#                             he_train=1,he_test=0,beta_sort_train=0,beta_sort_test=0,norm_type=norm_type,
#                             save_gene_data=F, save_surv_data=F, save_beta=F, plot_km=F, run_analysis=T,
#                             seed=1234)
# 
#   # Squared terms with intercept
#   nl.shiftqua.out = sim.survdata(surv.data,
#                                  N=10000,N.test=10000,test.ids=test.ids,sim.dataType='nl-shiftquad',
#                                  nonzero.genes=nonzero.genes,n_batch=10,h0=h0.shift,
#                                  he_train=1,he_test=0,beta_sort_train=0,beta_sort_test=0,norm_type=norm_type,
#                                  save_gene_data=F, save_surv_data=F, save_beta=F, plot_km=F, run_analysis=T,
#                                  seed=1234)
#   # Gene-gene interaction terms
#   nl.interact.out = sim.survdata(surv.data,
#                                  N=10000, N.test=10000, test.ids=test.ids, sim.dataType='nl-interaction',
#                                  nonzero.genes=nonzero.genes, n_batch=10, h0=h0.inter,
#                                  he_train=1,he_test=0,beta_sort_train=0,beta_sort_test=0,norm_type=norm_type,
#                                  save_gene_data=F, save_surv_data=F, save_beta=F, plot_km=F, run_analysis=T,
#                                  seed=1234)
#   # Sine interaction
#   nl.sine.out = sim.survdata(surv.data,
#                              N=10000, N.test=10000, test.ids=test.ids, sim.dataType='nl-sine',
#                              nonzero.genes=nonzero.genes, n_batch=10, h0=h0.sine,
#                              he_train=1,he_test=0,beta_sort_train=0,beta_sort_test=0,norm_type=norm_type,
#                              save_gene_data=F, save_surv_data=F, save_beta=F, plot_km=F, run_analysis=T,
#                              seed=1234)
# }
# 
# 
# 
## ################################# BE10Asso10 #####################################
for (norm_type in 0:6) {

  h0.mod <- 2.75
  h0.weak <- 65
  h0.qua <- 22
  h0.shift <- 19
  h0.inter <- 25
  h0.sine <- 0.8

  if (norm_type==5) {
    h0.weak = 9.5
    h0.qua = 15
    h0.shift = 16
  }

  # Linear risk with moderate effects
  linear.moderate.out = sim.survdata(surv.data,
                                     N=10000, N.test=10000, test.ids=test.ids, sim.dataType='linear-moderate',
                                     nonzero.genes=nonzero.genes, n_batch=10, h0=h0.mod,
                                     he_train=1,he_test=0,beta_sort_train=0.01,beta_sort_test=0,norm_type=norm_type,
                                     save_gene_data=F, save_surv_data=F, save_beta=F, plot_km=F, run_analysis=T,
                                     subset_sizes = c(100,500,1000,2000,5000,10000),
                                     runs_per_size = c(20,20,20,20,10,10),
                                     splits_per_size = c(3,5,5,10,10,10),
                                     seed=1234)

  # Linear risk with weak effects
  linear.weak.out = sim.survdata(surv.data,
                                 N=10000, N.test=10000, test.ids=test.ids, sim.dataType='linear-weak',
                                 nonzero.genes=nonzero.genes, n_batch=10, h0=h0.weak,
                                 he_train=1,he_test=0,beta_sort_train=0.01,beta_sort_test=0,norm_type=norm_type,
                                 save_gene_data=F, save_surv_data=F, save_beta=F, plot_km=F, run_analysis=T,
                                 subset_sizes = c(100,500,1000,2000,5000,10000),
                                 runs_per_size = c(20,20,20,20,10,10),
                                 splits_per_size = c(3,5,5,10,10,10),
                                 seed=1234)

  # Squared terms
  nl.qua.out = sim.survdata(surv.data,
                            N=10000, N.test=10000, test.ids=test.ids, sim.dataType='nl-quadratic',
                            nonzero.genes=nonzero.genes, n_batch=10, h0=h0.qua,
                            he_train=1,he_test=0,beta_sort_train=0.01,beta_sort_test=0,norm_type=norm_type,
                            save_gene_data=F, save_surv_data=F, save_beta=F, plot_km=F, run_analysis=T,
                            subset_sizes = c(100,500,1000,2000,5000,10000),
                            runs_per_size = c(20,20,20,20,10,10),
                            splits_per_size = c(3,5,5,10,10,10),
                            seed=1234)

  # Squared terms with intercept
  nl.shiftqua.out = sim.survdata(surv.data,
                                 N=10000,N.test=10000,test.ids=test.ids,sim.dataType='nl-shiftquad',
                                 nonzero.genes=nonzero.genes,n_batch=10,h0=h0.shift,
                                 he_train=1,he_test=0,beta_sort_train=0.01,beta_sort_test=0,norm_type=norm_type,
                                 save_gene_data=F, save_surv_data=F, save_beta=F, plot_km=F, run_analysis=T,
                                 subset_sizes = c(100,500,1000,2000,5000,10000),
                                 runs_per_size = c(20,20,20,20,10,10),
                                 splits_per_size = c(3,5,5,10,10,10),
                                 seed=1234)

  # Gene-gene interaction terms
  nl.interact.out = sim.survdata(surv.data,
                                 N=10000, N.test=10000, test.ids=test.ids, sim.dataType='nl-interaction',
                                 nonzero.genes=nonzero.genes, n_batch=10, h0=h0.inter,
                                 he_train=1,he_test=0,beta_sort_train=0.01,beta_sort_test=0,norm_type=norm_type,
                                 save_gene_data=F, save_surv_data=F, save_beta=F, plot_km=F, run_analysis=T,
                                 subset_sizes = c(100,500,1000,2000,5000,10000),
                                 runs_per_size = c(20,20,20,20,10,10),
                                 splits_per_size = c(3,5,5,10,10,10),
                                 seed=1234)

  # Sine interaction
  nl.sine.out = sim.survdata(surv.data,
                             N=10000, N.test=10000, test.ids=test.ids, sim.dataType='nl-sine',
                             nonzero.genes=nonzero.genes, n_batch=10, h0=h0.sine,
                             he_train=1,he_test=0,beta_sort_train=0.01,beta_sort_test=0,norm_type=norm_type,
                             save_gene_data=F, save_surv_data=F, save_beta=F, plot_km=F, run_analysis=T,
                             subset_sizes = c(100,500,1000,2000,5000,10000),
                             runs_per_size = c(20,20,20,20,10,10),
                             splits_per_size = c(3,5,5,10,10,10),
                             seed=1234)
}


# ## ################################# BE11Asso00 #####################################
# for (norm_type in 0:6) {
# 
#   h0.mod <- 2.75
#   h0.weak <- 65
#   h0.qua <- 22
#   h0.shift <- 19
#   h0.inter <- 25
#   h0.sine <- 0.8
# 
#   if (norm_type==5) {
#     h0.weak = 9.5
#     h0.qua = 15
#     h0.shift = 16
#   }
# 
#   # Linear risk with moderate effects
#   linear.moderate.out = sim.survdata(surv.data,
#                                      N=10000, N.test=10000, test.ids=test.ids, sim.dataType='linear-moderate',
#                                      nonzero.genes=nonzero.genes, n_batch=10, h0=h0.mod,
#                                      he_train=1,he_test=1,beta_sort_train=0,beta_sort_test=0,norm_type=norm_type,
#                                      save_gene_data=F, save_surv_data=F, save_beta=F, plot_km=F, run_analysis=T,
#                                      seed=1234)
# 
# 
#   # Linear risk with weak effects
#   linear.weak.out = sim.survdata(surv.data,
#                                  N=10000, N.test=10000, test.ids=test.ids, sim.dataType='linear-weak',
#                                  nonzero.genes=nonzero.genes, n_batch=10, h0=h0.weak,
#                                  he_train=1,he_test=1,beta_sort_train=0,beta_sort_test=0,norm_type=norm_type,
#                                  save_gene_data=F, save_surv_data=F, save_beta=F, plot_km=F, run_analysis=T,
# 
#                                  seed=1234)
# 
#   # Squared terms
#   nl.qua.out = sim.survdata(surv.data,
#                             N=10000, N.test=10000, test.ids=test.ids, sim.dataType='nl-quadratic',
#                             nonzero.genes=nonzero.genes, n_batch=10, h0=h0.qua,
#                             he_train=1,he_test=1,beta_sort_train=0,beta_sort_test=0,norm_type=norm_type,
#                             save_gene_data=F, save_surv_data=F, save_beta=F, plot_km=F, run_analysis=T,
#                             seed=1234)
# 
#   # Squared terms with intercept
#   nl.shiftqua.out = sim.survdata(surv.data,
#                                  N=10000,N.test=10000,test.ids=test.ids,sim.dataType='nl-shiftquad',
#                                  nonzero.genes=nonzero.genes,n_batch=10,h0=h0.shift,
#                                  he_train=1,he_test=1,beta_sort_train=0,beta_sort_test=0,norm_type=norm_type,
#                                  save_gene_data=F, save_surv_data=F, save_beta=F, plot_km=F, run_analysis=T,
#                                  seed=1234)
# 
#   # Gene-gene interaction terms
#   nl.interact.out = sim.survdata(surv.data,
#                                  N=10000, N.test=10000, test.ids=test.ids, sim.dataType='nl-interaction',
#                                  nonzero.genes=nonzero.genes, n_batch=10, h0=h0.inter,
#                                  he_train=1,he_test=1,beta_sort_train=0,beta_sort_test=0,norm_type=norm_type,
#                                  save_gene_data=F, save_surv_data=F, save_beta=F, plot_km=F, run_analysis=T,
#                                  seed=1234)
# 
#   # Sine interaction
#   nl.sine.out = sim.survdata(surv.data,
#                              N=10000, N.test=10000, test.ids=test.ids, sim.dataType='nl-sine',
#                              nonzero.genes=nonzero.genes, n_batch=10, h0=h0.sine,
#                              he_train=1,he_test=1,beta_sort_train=0,beta_sort_test=0,norm_type=norm_type,
#                              save_gene_data=F, save_surv_data=F, save_beta=F, plot_km=F, run_analysis=T,
#                              seed=1234)
# }




# ## ################################# BE11Asso10 #####################################
# for (norm_type in 0:6) {
#   
#   h0.mod <- 2.75
#   h0.weak <- 65
#   h0.qua <- 22
#   h0.shift <- 19
#   h0.inter <- 25
#   h0.sine <- 0.8
#   
#   if (norm_type==5) {
#     h0.weak = 9.5
#     h0.qua = 15
#     h0.shift = 16
#   }
#   
#   # Linear risk with moderate effects
#   linear.moderate.out = sim.survdata(surv.data,
#                                      N=10000, N.test=10000, test.ids=test.ids, sim.dataType='linear-moderate',
#                                      nonzero.genes=nonzero.genes, n_batch=10, h0=h0.mod,
#                                      he_train=1,he_test=1,beta_sort_train=0.05,beta_sort_test=0,norm_type=norm_type,
#                                      save_gene_data=F, save_surv_data=F, save_beta=F, plot_km=F, run_analysis=T,
#                                      seed=1234)
#   
#   # Linear risk with weak effects
#   linear.weak.out = sim.survdata(surv.data,
#                                  N=10000, N.test=10000, test.ids=test.ids, sim.dataType='linear-weak',
#                                  nonzero.genes=nonzero.genes, n_batch=10, h0=h0.weak,
#                                  he_train=1,he_test=1,beta_sort_train=0.05,beta_sort_test=0,norm_type=norm_type,
#                                  save_gene_data=F, save_surv_data=F, save_beta=F, plot_km=F, run_analysis=T,
#                                  seed=1234)
#   
#   # Squared terms
#   nl.qua.out = sim.survdata(surv.data,
#                             N=10000, N.test=10000, test.ids=test.ids, sim.dataType='nl-quadratic',
#                             nonzero.genes=nonzero.genes, n_batch=10, h0=h0.qua,
#                             he_train=1,he_test=1,beta_sort_train=0.05,beta_sort_test=0,norm_type=norm_type,
#                             save_gene_data=F, save_surv_data=F, save_beta=F, plot_km=F, run_analysis=T,
#                             seed=1234)
#   
#   # Squared terms with intercept
#   nl.shiftqua.out = sim.survdata(surv.data,
#                                  N=10000,N.test=10000,test.ids=test.ids,sim.dataType='nl-shiftquad',
#                                  nonzero.genes=nonzero.genes,n_batch=10,h0=h0.shift,
#                                  he_train=1,he_test=1,beta_sort_train=0.05,beta_sort_test=0,norm_type=norm_type,
#                                  save_gene_data=F, save_surv_data=F, save_beta=F, plot_km=F, run_analysis=T,
#                                  seed=1234)
#   
#   # Gene-gene interaction terms
#   nl.interact.out = sim.survdata(surv.data,
#                                  N=10000, N.test=10000, test.ids=test.ids, sim.dataType='nl-interaction',
#                                  nonzero.genes=nonzero.genes, n_batch=10, h0=h0.inter,
#                                  he_train=1,he_test=1,beta_sort_train=0.05,beta_sort_test=0,norm_type=norm_type,
#                                  save_gene_data=F, save_surv_data=F, save_beta=F, plot_km=F, run_analysis=T,
#                                  seed=1234)
#   
#   # Sine interaction
#   nl.sine.out = sim.survdata(surv.data,
#                              N=10000, N.test=10000, test.ids=test.ids, sim.dataType='nl-sine',
#                              nonzero.genes=nonzero.genes, n_batch=10, h0=h0.sine,
#                              he_train=1,he_test=1,beta_sort_train=0.05,beta_sort_test=0,norm_type=norm_type,
#                              save_gene_data=F, save_surv_data=F, save_beta=F, plot_km=F, run_analysis=T,
#                              seed=1234)
# }



# ## ################################# BE11Asso11 #####################################
# for (norm_type in 1:6) {
# 
#   h0.mod <- 2.75
#   h0.weak <- 65
#   h0.qua <- 22
#   h0.shift <- 19
#   h0.inter <- 25
#   h0.sine <- 0.8
# 
#   if (norm_type==5) {
#     h0.weak = 9.5
#     h0.qua = 15
#     h0.shift = 16
#   }
# 
#   # Linear risk with moderate effects
#   linear.moderate.out = sim.survdata(surv.data,
#                                      N=10000, N.test=10000, test.ids=test.ids, sim.dataType='linear-moderate',
#                                      nonzero.genes=nonzero.genes, n_batch=10, h0=h0.mod,
#                                      he_train=1,he_test=1,beta_sort_train=0.05,beta_sort_test=0.05,norm_type=norm_type,
#                                      save_gene_data=F, save_surv_data=F, save_beta=F, plot_km=F, run_analysis=T,
#                                      seed=1234)
# 
#   # Linear risk with weak effects
#   linear.weak.out = sim.survdata(surv.data,
#                                  N=10000, N.test=10000, test.ids=test.ids, sim.dataType='linear-weak',
#                                  nonzero.genes=nonzero.genes, n_batch=10, h0=h0.weak,
#                                  he_train=1,he_test=1,beta_sort_train=0.05,beta_sort_test=0.05,norm_type=norm_type,
#                                  save_gene_data=F, save_surv_data=F, save_beta=F, plot_km=F, run_analysis=T,
#                                  seed=1234)
# 
#   # Squared terms
#   nl.qua.out = sim.survdata(surv.data,
#                             N=10000, N.test=10000, test.ids=test.ids, sim.dataType='nl-quadratic',
#                             nonzero.genes=nonzero.genes, n_batch=10, h0=h0.qua,
#                             he_train=1,he_test=1,beta_sort_train=0.05,beta_sort_test=0.05,norm_type=norm_type,
#                             save_gene_data=F, save_surv_data=F, save_beta=F, plot_km=F, run_analysis=T,
#                             seed=1234)
# 
#   # Squared terms with intercept
#   nl.shiftqua.out = sim.survdata(surv.data,
#                                  N=10000,N.test=10000,test.ids=test.ids,sim.dataType='nl-shiftquad',
#                                  nonzero.genes=nonzero.genes,n_batch=10,h0=h0.shift,
#                                  he_train=1,he_test=1,beta_sort_train=0.05,beta_sort_test=0.05,norm_type=norm_type,
#                                  save_gene_data=F, save_surv_data=F, save_beta=F, plot_km=F, run_analysis=T,
#                                  seed=1234)
# 
#   # Gene-gene interaction terms
#   nl.interact.out = sim.survdata(surv.data,
#                                  N=10000, N.test=10000, test.ids=test.ids, sim.dataType='nl-interaction',
#                                  nonzero.genes=nonzero.genes, n_batch=10, h0=h0.inter,
#                                  he_train=1,he_test=1,beta_sort_train=0.05,beta_sort_test=0.05,norm_type=norm_type,
#                                  save_gene_data=F, save_surv_data=F, save_beta=F, plot_km=F, run_analysis=T,
#                                  seed=1234)
# 
#   # Sine interaction
#   nl.sine.out = sim.survdata(surv.data,
#                              N=10000, N.test=10000, test.ids=test.ids, sim.dataType='nl-sine',
#                              nonzero.genes=nonzero.genes, n_batch=10, h0=h0.sine,
#                              he_train=1,he_test=1,beta_sort_train=0.05,beta_sort_test=0.05,norm_type=norm_type,
#                              save_gene_data=F, save_surv_data=F, save_beta=F, plot_km=F, run_analysis=T,
#                              seed=1234)
# }
