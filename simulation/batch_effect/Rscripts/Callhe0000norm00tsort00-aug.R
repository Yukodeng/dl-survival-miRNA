source('~/dl-survival-miRNA/simulation/batch_effect/with_without_batch_maf/aug-mainSimSeqedgeRparV2.R')

seed <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))

mainSimSeqedgeRpar.aug(c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9),c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444),300,200,0,0,30,20,10,0,0,0,0,5,5,1,1,0,0,0,0,0,0,0,0.9,c(0,0.05),20,0.3,30,5,'he0000bs20norm00tsort00-aug',20,seed)
