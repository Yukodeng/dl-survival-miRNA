# options(echo=TRUE)

# NOTE: supply job array number 1-7
# norm_type <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID')) - 1

source('~/dl-survival-miRNA/simulation/aug-batch-survival-pipeline-v2.R')

# Run Simulation ---------------------------------

## ############################# BE00Asso00 ###############################
# for (i in 1:7) {
#   norm_type = i-1
#   # Linear risk with moderate effects
#   linear.moderate.out = sim.survdata(surv.data,
#                                      N=10000, N.test=10000, test.ids=test.ids, sim.dataType='linear-moderate',
#                                      nonzero.genes=nonzero.genes, n_batch=10, h0=2.75,
#                                      he_train=0,he_test=0,beta_sort_train=0,beta_sort_test=0,norm_type=norm_type,
#                                      save_gene_data=F, save_surv_data=F, save_beta=F, plot_km=F, run_analysis=F,
#                                      seed=1234)
#   # Linear risk with weak effects
#   linear.weak.out = sim.survdata(surv.data,
#                                  N=10000, N.test=10000, test.ids=test.ids, sim.dataType='linear-weak',
#                                  nonzero.genes=nonzero.genes, n_batch=10, h0=65,
#                                  he_train=0,he_test=0,beta_sort_train=0,beta_sort_test=0,norm_type=norm_type,
#                                  save_gene_data=F, save_surv_data=F, save_beta=F, plot_km=F, run_analysis=F,
#                                  seed=1234)
#   # Squared terms
#   nl.qua.out = sim.survdata(surv.data,
#                             N=10000, N.test=10000, test.ids=test.ids, sim.dataType='nl-quadratic',
#                             nonzero.genes=nonzero.genes, n_batch=10, h0=22,
#                             he_train=0,he_test=0,beta_sort_train=0,beta_sort_test=0,norm_type=norm_type,
#                             save_gene_data=F, save_surv_data=F, save_beta=F, plot_km=F, run_analysis=F,
#                             seed=1234)
#   # Squared terms with intercept
#   nl.shiftqua.out = sim.survdata(surv.data,
#                                  N=10000,N.test=10000,test.ids=test.ids,sim.dataType='nl-shiftquad',
#                                  nonzero.genes=nonzero.genes,n_batch=10,h0=19,
#                                  he_train=0,he_test=0,beta_sort_train=0,beta_sort_test=0,norm_type=norm_type,
#                                  save_gene_data=F,save_surv_data=F,save_beta=F,plot_km=F,run_analysis=F,
#                                  seed=1234)
#   # Gene-gene interaction terms
#   nl.interact.out = sim.survdata(surv.data,
#                                  N=10000, N.test=10000, test.ids=test.ids, sim.dataType='nl-interaction',
#                                  nonzero.genes=nonzero.genes, n_batch=10, h0=25,
#                                  he_train=0,he_test=0,beta_sort_train=0,beta_sort_test=0,norm_type=norm_type,
#                                  save_gene_data=F, save_surv_data=F, save_beta=F, plot_km=F, run_analysis=F,
#                                  seed=1234)
#   # Sine interaction
#   nl.sine.out = sim.survdata(surv.data,
#                              N=10000, N.test=10000, test.ids=test.ids, sim.dataType='nl-sine',
#                              nonzero.genes=nonzero.genes, n_batch=10, h0=.8,
#                              he_train=0,he_test=0,beta_sort_train=0,beta_sort_test=0,norm_type=norm_type,
#                              save_gene_data=F, save_surv_data=F, save_beta=F, plot_km=F, run_analysis=F,
#                              seed=1234)
# }

## ############################# BE10Asso00 ####################################
for (i in 1:7) {
  norm_type = i-1
  # Linear risk with moderate effects
  linear.moderate.out = sim.survdata(surv.data,
                                     N=10000, N.test=10000, test.ids=test.ids, sim.dataType='linear-moderate',
                                     nonzero.genes=nonzero.genes, n_batch=10, h0=2.75,
                                     he_train=1,he_test=0,beta_sort_train=0,beta_sort_test=0,norm_type=norm_type,
                                     save_gene_data=T, save_surv_data=T, save_beta=T, plot_km=F, run_analysis=T,
                                     seed=1234)
  # Linear risk with weak effects
  linear.weak.out = sim.survdata(surv.data,
                                 N=10000, N.test=10000, test.ids=test.ids, sim.dataType='linear-weak',
                                 nonzero.genes=nonzero.genes, n_batch=10, h0=65,
                                 he_train=1,he_test=0,beta_sort_train=0,beta_sort_test=0,norm_type=norm_type,
                                 save_gene_data=T, save_surv_data=T, save_beta=T, plot_km=F, run_analysis=T,
                                 seed=1234)
  # Squared terms
  nl.qua.out = sim.survdata(surv.data,
                            N=10000, N.test=10000, test.ids=test.ids, sim.dataType='nl-quadratic',
                            nonzero.genes=nonzero.genes, n_batch=10, h0=22,
                            he_train=1,he_test=0,beta_sort_train=0,beta_sort_test=0,norm_type=norm_type,
                            save_gene_data=T, save_surv_data=T, save_beta=T, plot_km=F, run_analysis=T,
                            seed=1234)
  # Squared terms with intercept
  nl.shiftqua.out = sim.survdata(surv.data,
                                 N=10000,N.test=10000,test.ids=test.ids,sim.dataType='nl-shiftquad',
                                 nonzero.genes=nonzero.genes,n_batch=10,h0=19,
                                 he_train=1,he_test=0,beta_sort_train=0,beta_sort_test=0,norm_type=norm_type,
                                 save_gene_data=T, save_surv_data=T, save_beta=T, plot_km=F, run_analysis=T,
                                 seed=1234)
  # Gene-gene interaction terms
  nl.interact.out = sim.survdata(surv.data,
                                 N=10000, N.test=10000, test.ids=test.ids, sim.dataType='nl-interaction',
                                 nonzero.genes=nonzero.genes, n_batch=10, h0=25,
                                 he_train=1,he_test=0,beta_sort_train=0,beta_sort_test=0,norm_type=norm_type,
                                 save_gene_data=T, save_surv_data=T, save_beta=T, plot_km=F, run_analysis=T,
                                 seed=1234)
  # Sine interaction
  nl.sine.out = sim.survdata(surv.data,
                             N=10000, N.test=10000, test.ids=test.ids, sim.dataType='nl-sine',
                             nonzero.genes=nonzero.genes, n_batch=10, h0=.8,
                             he_train=1,he_test=0,beta_sort_train=0,beta_sort_test=0,norm_type=norm_type,
                             save_gene_data=T, save_surv_data=T, save_beta=T, plot_km=F, run_analysis=T,
                             seed=1234)
}

## ################################# BE10Asso10 #####################################
# for (i in 1:7) {
#   norm_type = i-1
#   # Linear risk with moderate effects
#   linear.moderate.out = sim.survdata(surv.data,
#                                      N=10000, N.test=10000, test.ids=test.ids, sim.dataType='linear-moderate',
#                                      nonzero.genes=nonzero.genes, n_batch=10, h0=2.75,
#                                      he_train=1,he_test=0,beta_sort_train=0.05,beta_sort_test=0,norm_type=norm_type,
#                                      save_gene_data=T, save_surv_data=T, save_beta=T, plot_km=F, run_analysis=T,
#                                      seed=1234)
#   
#   # Linear risk with weak effects
#   linear.weak.out = sim.survdata(surv.data,
#                                  N=10000, N.test=10000, test.ids=test.ids, sim.dataType='linear-weak',
#                                  nonzero.genes=nonzero.genes, n_batch=10, h0=65,
#                                  he_train=1,he_test=0,beta_sort_train=0.05,beta_sort_test=0,norm_type=norm_type,
#                                  save_gene_data=F, save_surv_data=T, save_beta=T, plot_km=F, run_analysis=T,
#                                  seed=1234)
#   
#   # Squared terms
#   nl.qua.out = sim.survdata(surv.data,
#                             N=10000, N.test=10000, test.ids=test.ids, sim.dataType='nl-quadratic',
#                             nonzero.genes=nonzero.genes, n_batch=10, h0=22,
#                             he_train=1,he_test=0,beta_sort_train=0.05,beta_sort_test=0,norm_type=norm_type,
#                             save_gene_data=F, save_surv_data=T, save_beta=T, plot_km=F, run_analysis=T,
#                             seed=1234)
#   
#   # Squared terms with intercept
#   nl.shiftqua.out = sim.survdata(surv.data,
#                                  N=10000,N.test=10000,test.ids=test.ids,sim.dataType='nl-shiftquad',
#                                  nonzero.genes=nonzero.genes,n_batch=10,h0=19,
#                                  he_train=1,he_test=0,beta_sort_train=0.05,beta_sort_test=0,norm_type=norm_type,
#                                  save_gene_data=F, save_surv_data=T, save_beta=T, plot_km=F, run_analysis=T,
#                                  seed=1234)
#   
#   # Gene-gene interaction terms
#   nl.interact.out = sim.survdata(surv.data,
#                                  N=10000, N.test=10000, test.ids=test.ids, sim.dataType='nl-interaction',
#                                  nonzero.genes=nonzero.genes, n_batch=10, h0=25,
#                                  he_train=1,he_test=0,beta_sort_train=0.05,beta_sort_test=0,norm_type=norm_type,
#                                  save_gene_data=F, save_surv_data=T, save_beta=T, plot_km=F, run_analysis=T,
#                                  seed=1234)
#   
#   # Sine interaction
#   nl.sine.out = sim.survdata(surv.data,
#                              N=10000, N.test=10000, test.ids=test.ids, sim.dataType='nl-sine',
#                              nonzero.genes=nonzero.genes, n_batch=10, h0=.8,
#                              he_train=1,he_test=0,beta_sort_train=0.05,beta_sort_test=0,norm_type=norm_type,
#                              save_gene_data=F, save_surv_data=T, save_beta=T, plot_km=F, run_analysis=T,
#                              seed=1234)
# }

## ################################# BE11Asso00 #####################################
# for (i in 1:7) {
#   norm_type = i-1
#   # Linear risk with moderate effects
#   linear.moderate.out = sim.survdata(surv.data,
#                                      N=10000, N.test=10000, test.ids=test.ids, sim.dataType='linear-moderate',
#                                      nonzero.genes=nonzero.genes, n_batch=10, h0=2.75,
#                                      he_train=1,he_test=1,beta_sort_train=0,beta_sort_test=0,norm_type=norm_type,
#                                      save_gene_data=T, save_surv_data=T, save_beta=T, plot_km=F, run_analysis=T,
#                                      seed=1234)
#   
#   # Linear risk with weak effects
#   linear.weak.out = sim.survdata(surv.data,
#                                  N=10000, N.test=10000, test.ids=test.ids, sim.dataType='linear-weak',
#                                  nonzero.genes=nonzero.genes, n_batch=10, h0=65,
#                                  he_train=1,he_test=1,beta_sort_train=0,beta_sort_test=0,norm_type=norm_type,
#                                  save_gene_data=F, save_surv_data=T, save_beta=T, plot_km=F, run_analysis=T,
#                                  seed=1234)
#   
#   # Squared terms
#   nl.qua.out = sim.survdata(surv.data,
#                             N=10000, N.test=10000, test.ids=test.ids, sim.dataType='nl-quadratic',
#                             nonzero.genes=nonzero.genes, n_batch=10, h0=22,
#                             he_train=1,he_test=1,beta_sort_train=0,beta_sort_test=0,norm_type=norm_type,
#                             save_gene_data=F, save_surv_data=T, save_beta=T, plot_km=F, run_analysis=T,
#                             seed=1234)
#   
#   # Squared terms with intercept
#   nl.shiftqua.out = sim.survdata(surv.data,
#                                  N=10000,N.test=10000,test.ids=test.ids,sim.dataType='nl-shiftquad',
#                                  nonzero.genes=nonzero.genes,n_batch=10,h0=19,
#                                  he_train=1,he_test=1,beta_sort_train=0,beta_sort_test=0,norm_type=norm_type,
#                                  save_gene_data=F, save_surv_data=T, save_beta=T, plot_km=F, run_analysis=T,
#                                  seed=1234)
#   
#   # Gene-gene interaction terms
#   nl.interact.out = sim.survdata(surv.data,
#                                  N=10000, N.test=10000, test.ids=test.ids, sim.dataType='nl-interaction',
#                                  nonzero.genes=nonzero.genes, n_batch=10, h0=25,
#                                  he_train=1,he_test=1,beta_sort_train=0,beta_sort_test=0,norm_type=norm_type,
#                                  save_gene_data=F, save_surv_data=T, save_beta=T, plot_km=F, run_analysis=T,
#                                  seed=1234)
#   
#   # Sine interaction
#   nl.sine.out = sim.survdata(surv.data,
#                              N=10000, N.test=10000, test.ids=test.ids, sim.dataType='nl-sine',
#                              nonzero.genes=nonzero.genes, n_batch=10, h0=.8,
#                              he_train=1,he_test=1,beta_sort_train=0,beta_sort_test=0,norm_type=norm_type,
#                              save_gene_data=F, save_surv_data=T, save_beta=T, plot_km=F, run_analysis=T,
#                              seed=1234)
# }