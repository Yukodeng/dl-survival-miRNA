resultproc=function(nonzero_position, b0value, p, stratify, cmbt, nbatch, exclude, addon, output_name){
  #### [1/27/2025] (YD) added the following 3 lines
  output_path=paste0("~/dl-survival-miRNA/simulation/batch_effect/data/",addon)
  setwd(output_path)
  dir.create(file.path(output_path, 'organized'), showWarnings = F)
  
  b0=rep(0,p)
  b0[nonzero_position]=b0value
  
  for(second in c("o","u","l","a")){
    bname=paste("b_",second,sep=""); btemp=c()
    cname_test=paste("c_",second,"_test",sep=""); ctemp_test=c()
    if(second %in% c("l","a")){
      lname=paste("l_",second,sep=""); ltemp=c()
    }      
    # pssizes=c()
    for(i in (1:nbatch)[-exclude]){
      btempi=read.csv(paste(bname, "_", addon, "_", i, ".csv", sep=""))
      btemp=rbind(btemp, btempi)
      ctempi_test=read.csv(paste(cname_test, "_", addon, "_", i, ".csv", sep=""))
      ctemp_test=c(ctemp_test, ctempi_test$x)
      if(second %in% c("l","a")){
        ltempi=read.csv(paste(lname, "_", addon, "_", i, ".csv", sep=""))
        ltemp=c(ltemp, ltempi$x)
      }
      # pssizei=read.csv(paste("pssize_", addon, "_", i, ".csv", sep=""))
      # pssizes=c(pssizes, pssizei)
    }
    assign(bname, btemp)
    assign(cname_test, ctemp_test)
    if(second %in% c("l","a")){
      assign(lname, ltemp)
    }
  }

  p0=sum(b0==0)
  p1=sum(b0!=0)
  b0mat=rep(1,nrow(btemp))%*%t(b0)
  
  cs_test=cs_ciu_test=cs_cil_test=c() # mean Harrell's C statistic from test data
  l2error=l2error_ciu=l2error_cil=c() # mean estimation error of betas (mean squared error loss, i.e. L2 norm)
  ze0=c() # percent of correctly and incorrectly identified 0 betas
  ze=c()  # percent of correctly and incorrectly identified nonzero betas
  pc=c()  # percent of correctly identified parameter
  fdr=c() # false discovery rate
  ritm100=ritm99=c() # different rates of identifying the true model (RITM)
  size=size_ciu=size_cil=c() # median model size
  
  for(second in c("o","u","l","a")){
    bname=paste("b_",second,sep="")
    cname_test=paste("c_",second,"_test",sep="")
    bfull=eval(parse(text=bname))
    cfull_test=eval(parse(text=cname_test))
    cs_test=c(cs_test, mean(cfull_test, na.rm=T))
    cs_ciu_test=c(cs_ciu_test, quantile(cfull_test, 0.975, na.rm=T))
    cs_cil_test=c(cs_cil_test, quantile(cfull_test, 0.025, na.rm=T))
    l2error=c(l2error, mean(sqrt(diag(t((t(bfull)-b0))%*%(t(bfull)-b0))), na.rm=T))
    l2error_ciu=c(l2error_ciu, quantile(sqrt(diag(t((t(bfull)-b0))%*%(t(bfull)-b0))), 0.975, na.rm=T))
    l2error_cil=c(l2error_cil, quantile(sqrt(diag(t((t(bfull)-b0))%*%(t(bfull)-b0))), 0.025, na.rm=T))
    ze0=rbind(ze0, c(mean(rowSums(bfull[,b0==0]==0, na.rm=T))/p0, mean(rowSums(bfull[,b0!=0]==0, na.rm=T))/p1))
    ze=rbind(ze, c(mean(rowSums(bfull[,b0!=0]!=0, na.rm=T))/p1, mean(rowSums(bfull[,b0==0]!=0, na.rm=T))/p0))
    pc=c(pc, mean(rowSums(bfull[,b0==0]==0, na.rm=T)+rowSums(bfull[,b0!=0]!=0, na.rm=T))/p)
    fdr=c(fdr, mean(rowSums(b0mat==0 & bfull!=0, na.rm=T)/rowSums(bfull!=0, na.rm=T), na.rm=T))
    ritm100=c(ritm100, mean((rowSums(bfull[,b0==0]==0, na.rm=T)+rowSums(bfull[,b0!=0]!=0, na.rm=T))==p)*100)
    ritm99=c(ritm99, mean((rowSums(bfull[,b0==0]==0, na.rm=T)+rowSums(bfull[,b0!=0]!=0, na.rm=T))>=(0.99*p))*100)
    size=c(size, median(rowSums(bfull!=0, na.rm=T)))
    size_ciu=c(size_ciu, quantile(rowSums(bfull!=0), 0.975, na.rm=T))
    size_cil=c(size_cil, quantile(rowSums(bfull!=0), 0.025, na.rm=T))
  }

  master_result=cbind(cs_test, cs_cil_test, cs_ciu_test, l2error, l2error_cil, l2error_ciu, 
                      ze0, ze, pc, fdr, ritm100, ritm99, size, size_cil, size_ciu)

  rownames(master_result)=c("o","u","l","a")
  colnames(master_result)=c("cs_test","cs_cil_test","cs_ciu_test","l2error","l2error_cil","l2error_ciu",
                            "TN","FN","TP","FP","pc","FDR","ritm100","ritm99","size","size_cil","size_ciu")
  write.csv(master_result, file.path('organized', paste(output_name, "_results.csv", sep="")))
  
  # write.csv(pssizes, paste("pssize_",addon,".csv",sep=""), row.names=F)
  
  for(second in c("o","u","l","a")){
    bname=paste("b_",second,sep="")
    bfull=eval(parse(text=bname))
    write.csv(bfull, file.path('organized', paste(bname,"_",addon,".csv",sep="")), row.names=F)
    
    cname_test=paste("c_",second,"_test",sep="")
    cfull_test=eval(parse(text=cname_test))
    write.csv(cfull_test, file.path('organized',paste(cname_test,"_",addon,".csv",sep="")), row.names=F)
    
    if(second %in% c("l","a")){
      lname=paste("l_",second,sep="")
      lfull=eval(parse(text=lname))
      write.csv(lfull, file.path('organized', paste(lname,"_",addon,".csv",sep="")), row.names=F)
    }
  }

  
  if(stratify==1){
    for(second in c("o","u","l","a")){
      bname=paste("b_",second,"s",sep=""); btemp=c()
      cname_test=paste("c_",second,"s_test",sep=""); ctemp_test=c()
      if(second %in% c("l","a")){
        lname=paste("l_",second,"s",sep=""); ltemp=c()
      }
      for(i in (1:nbatch)[-exclude]){
        btempi=read.csv(paste(bname, "_", addon, "_", i, ".csv", sep=""))
        btemp=rbind(btemp, btempi)
        ctempi_test=read.csv(paste(cname_test, "_", addon, "_", i, ".csv", sep=""))
        ctemp_test=c(ctemp_test, ctempi_test$x)
        if(second %in% c("l","a")){
          ltempi=read.csv(paste(lname, "_", addon, "_", i, ".csv", sep=""))
          ltemp=c(ltemp, ltempi$x)  
        }
      }
      assign(bname, btemp)
      assign(cname_test, ctemp_test)
      if(second %in% c("l","a")){
        assign(lname, ltemp)
      }
    }
    
    p0=sum(b0==0)
    p1=sum(b0!=0)
    
    cs_test=cs_ciu_test=cs_cil_test=c() # mean Harrell's C statistic from test data
    l2error=l2error_ciu=l2error_cil=c() # mean estimation error of betas (mean squared error loss, i.e. L2 norm)
    ze0=c() # percent of correctly and incorrectly identified 0 betas
    ze=c()  # percent of correctly and incorrectly identified nonzero betas
    pc=c()  # percent of correctly identified parameter
    fdr=c() # false discovery rate
    ritm100=ritm99=c() # different rates of identifying the true model (RITM)
    size=size_ciu=size_cil=c() # median model size
    
    for(second in c("o","u","l","a")){
      bname=paste("b_",second,"s",sep="")
      cname_test=paste("c_",second,"s_test",sep="")
      bfull=eval(parse(text=bname))
      cfull_test=eval(parse(text=cname_test))
      cs_test=c(cs_test, mean(cfull_test, na.rm=T))
      cs_ciu_test=c(cs_ciu_test, quantile(cfull_test, 0.975, na.rm=T))
      cs_cil_test=c(cs_cil_test, quantile(cfull_test, 0.025, na.rm=T))
      l2error=c(l2error, mean(sqrt(diag(t((t(bfull)-b0))%*%(t(bfull)-b0))), na.rm=T))
      l2error_ciu=c(l2error_ciu, quantile(sqrt(diag(t((t(bfull)-b0))%*%(t(bfull)-b0))), 0.975, na.rm=T))
      l2error_cil=c(l2error_cil, quantile(sqrt(diag(t((t(bfull)-b0))%*%(t(bfull)-b0))), 0.025, na.rm=T))
      ze0=rbind(ze0, c(mean(rowSums(bfull[,b0==0]==0, na.rm=T))/p0, mean(rowSums(bfull[,b0!=0]==0, na.rm=T))/p1))
      ze=rbind(ze, c(mean(rowSums(bfull[,b0!=0]!=0, na.rm=T))/p1, mean(rowSums(bfull[,b0==0]!=0, na.rm=T))/p0))
      pc=c(pc, mean(rowSums(bfull[,b0==0]==0, na.rm=T)+rowSums(bfull[,b0!=0]!=0, na.rm=T))/p)
      fdr=c(fdr, mean(rowSums(b0mat==0 & bfull!=0, na.rm=T)/rowSums(bfull!=0, na.rm=T), na.rm=T))
      ritm100=c(ritm100, mean((rowSums(bfull[,b0==0]==0, na.rm=T)+rowSums(bfull[,b0!=0]!=0, na.rm=T))==p)*100)
      ritm99=c(ritm99, mean((rowSums(bfull[,b0==0]==0, na.rm=T)+rowSums(bfull[,b0!=0]!=0, na.rm=T))>=(0.99*p))*100)
      size=c(size, median(rowSums(bfull!=0, na.rm=T)))
      size_ciu=c(size_ciu, quantile(rowSums(bfull!=0), 0.975, na.rm=T))
      size_cil=c(size_cil, quantile(rowSums(bfull!=0), 0.025, na.rm=T))
    }
    
    master_result=cbind(cs_test, cs_cil_test, cs_ciu_test, l2error, l2error_cil, l2error_ciu, 
                        ze0, ze, pc, fdr, ritm100, ritm99, size, size_cil, size_ciu)
    rownames(master_result)=c("o","u","l","a")
    colnames(master_result)=c("cs_test","cs_cil_test","cs_ciu_test","l2error","l2error_cil","l2error_ciu",
                              "TN","FN","TP","FP","pc","FDR","ritm100","ritm99","size","size_cil","size_ciu")
    write.csv(master_result, file.path('organized', paste(output_name, "_strat_results.csv", sep="")))
    
    for(second in c("o","u","l","a")){
      bname=paste("b_",second,"s",sep="")
      bfull=eval(parse(text=bname))
      write.csv(bfull, file.path('organized', paste(bname,"_",addon,".csv",sep="")), row.names=F)
      
      cname_test=paste("c_",second,"s_test",sep="")
      cfull_test=eval(parse(text=cname_test))
      write.csv(cfull_test, file.path('organized', paste(cname_test,"_",addon,".csv",sep="")), row.names=F)
      
      if(second %in% c("l","a")){
        lname=paste("l_",second,"s",sep="")
        lfull=eval(parse(text=lname))
        write.csv(lfull, file.path('organized',paste(lname,"_",addon,".csv",sep="")), row.names=F)  
      }
    }
  }
  
  
  if(cmbt==1){
    for(second in c("o","u","l","a")){
      bname=paste("b_",second,"c",sep=""); btemp=c()
      cname_test=paste("c_",second,"c_test",sep=""); ctemp_test=c()
      if(second %in% c("l","a")){
        lname=paste("l_",second,"c",sep=""); ltemp=c()
      }
      for(i in (1:nbatch)[-exclude]){
        btempi=read.csv(paste(bname, "_", addon, "_", i, ".csv", sep=""))
        btemp=rbind(btemp, btempi)
        ctempi_test=read.csv(paste(cname_test, "_", addon, "_", i, ".csv", sep=""))
        ctemp_test=c(ctemp_test, ctempi_test$x)
        if(second %in% c("l","a")){
          ltempi=read.csv(paste(lname, "_", addon, "_", i, ".csv", sep=""))
          ltemp=c(ltemp, ltempi$x)  
        }
      }
      assign(bname, btemp)
      assign(cname_test, ctemp_test)
      if(second %in% c("l","a")){
        assign(lname, ltemp)
      }
    }
    
    p0=sum(b0==0)
    p1=sum(b0!=0)
    
    cs_test=cs_ciu_test=cs_cil_test=c() # mean Harrell's C statistic from test data
    l2error=l2error_ciu=l2error_cil=c() # mean estimation error of betas (mean squared error loss, i.e. L2 norm)
    ze0=c() # percent of correctly and incorrectly identified 0 betas
    ze=c()  # percent of correctly and incorrectly identified nonzero betas
    pc=c()  # percent of correctly identified parameter
    fdr=c() # false discovery rate
    ritm100=ritm99=c() # different rates of identifying the true model (RITM)
    size=size_ciu=size_cil=c() # median model size
    
    for(second in c("o","u","l","a")){
      bname=paste("b_",second,"c",sep="")
      cname_test=paste("c_",second,"c_test",sep="")
      bfull=eval(parse(text=bname))
      cfull_test=eval(parse(text=cname_test))
      cs_test=c(cs_test, mean(cfull_test, na.rm=T))
      cs_ciu_test=c(cs_ciu_test, quantile(cfull_test, 0.975, na.rm=T))
      cs_cil_test=c(cs_cil_test, quantile(cfull_test, 0.025, na.rm=T))
      l2error=c(l2error, mean(sqrt(diag(t((t(bfull)-b0))%*%(t(bfull)-b0))), na.rm=T))
      l2error_ciu=c(l2error_ciu, quantile(sqrt(diag(t((t(bfull)-b0))%*%(t(bfull)-b0))), 0.975, na.rm=T))
      l2error_cil=c(l2error_cil, quantile(sqrt(diag(t((t(bfull)-b0))%*%(t(bfull)-b0))), 0.025, na.rm=T))
      ze0=rbind(ze0, c(mean(rowSums(bfull[,b0==0]==0, na.rm=T))/p0, mean(rowSums(bfull[,b0!=0]==0, na.rm=T))/p1))
      ze=rbind(ze, c(mean(rowSums(bfull[,b0!=0]!=0, na.rm=T))/p1, mean(rowSums(bfull[,b0==0]!=0, na.rm=T))/p0))
      pc=c(pc, mean(rowSums(bfull[,b0==0]==0, na.rm=T)+rowSums(bfull[,b0!=0]!=0, na.rm=T))/p)
      fdr=c(fdr, mean(rowSums(b0mat==0 & bfull!=0, na.rm=T)/rowSums(bfull!=0, na.rm=T), na.rm=T))
      ritm100=c(ritm100, mean((rowSums(bfull[,b0==0]==0, na.rm=T)+rowSums(bfull[,b0!=0]!=0, na.rm=T))==p)*100)
      ritm99=c(ritm99, mean((rowSums(bfull[,b0==0]==0, na.rm=T)+rowSums(bfull[,b0!=0]!=0, na.rm=T))>=(0.99*p))*100)
      size=c(size, median(rowSums(bfull!=0, na.rm=T)))
      size_ciu=c(size_ciu, quantile(rowSums(bfull!=0), 0.975, na.rm=T))
      size_cil=c(size_cil, quantile(rowSums(bfull!=0), 0.025, na.rm=T))
    }
    
    master_result=cbind(cs_test, cs_cil_test, cs_ciu_test, l2error, l2error_cil, l2error_ciu, 
                         ze0, ze, pc, fdr, ritm100, ritm99, size, size_cil, size_ciu)
    rownames(master_result)=c("o","u","l","a")
    colnames(master_result)=c("cs_test","cs_cil_test","cs_ciu_test","l2error","l2error_cil","l2error_ciu",
                              "TN","FN","TP","FP","pc","FDR","ritm100","ritm99","size","size_cil","size_ciu")
    write.csv(master_result, file.path('organized', paste(output_name, "_combat_results.csv", sep="")))
    
    for(second in c("o","u","l","a")){
      bname=paste("b_",second,"c",sep="")
      bfull=eval(parse(text=bname))
      write.csv(bfull, file.path('organized', paste(bname,"_",addon,".csv",sep="")), row.names=F)
      
      cname_test=paste("c_",second,"c_test",sep="")
      cfull_test=eval(parse(text=cname_test))
      write.csv(cfull_test, file.path('organized', paste(cname_test,"_",addon,".csv",sep="")), row.names=F)
      
      if(second %in% c("l","a")){
        lname=paste("l_",second,"c",sep="")
        lfull=eval(parse(text=lname))
        write.csv(lfull, file.path('organized', paste(lname,"_",addon,".csv",sep="")), row.names=F)  
      }
    }
  }
}









