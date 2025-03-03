source('~/dl-survival-miRNA/simulation/batch_effect/mainSimSeqedgeRparV2.R')

seed <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))

mainSimSeqedgeRpar(c(2.4, -2.6, 2.0, -2.3, 5.0, -2.3, 1.8, -4.3, 3.1, -1.9),c(50, 51, 67, 85, 97, 123, 162, 260, 401, 444),300,200,0.9,0,0,20,0,0,0,0,5,5,1,1,0,0,0,2,1,1,0,c(0,0.05),20,0.3,30,5,'he0000bs20norm22tsort00',20, seed)
