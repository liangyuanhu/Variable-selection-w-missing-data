perform_met<-function(y, yhat, phat){
  y_new<-ifelse(is.na(y)==0 & y=="No",0,
                ifelse(is.na(y)==0 & y=="Yes",1,NA))
  accu<-confusionMatrix(yhat, y)$overall["Accuracy"]
  kappa<-confusionMatrix(yhat, y)$overall["Kappa"]
  sen<-confusionMatrix(yhat, y)$byClass["Sensitivity"]
  spe<-confusionMatrix(yhat, y)$byClass["Specificity"]
  AUC<-auc(roc(response=y_new, predictor=phat, direction="<"))
  perform<-matrix(c(AUC, sen, spe, accu, kappa),nrow=1)
  colnames(perform)<-c("auc","sen","spe","accu","kapppa")
  return(perform)
}
auc_met<-function(y, yhat, phat){
  y_new<-ifelse(is.na(y)==0 & y=="No",0,
                ifelse(is.na(y)==0 & y=="Yes",1,NA))
  AUC<-auc(roc(response=y_new, predictor=phat, direction="<"))
  return(AUC)
}
class_error<-function(y, yhat, phat=NULL){
  Metrics::ce(actual=y, predicted=yhat)
}


steplogit_p<-function(datax, event, direction, alpha_in=0.05, alpha_out=0.05,
                      exclude=NULL){
  # Exclude variables not use in the selection
  varname<-names(datax)[-which(names(datax)%in%c(event, exclude))]
  data0<-datax[,c(event,varname)]
  
  #
  if(direction%in%c("for_back","forward")){
    FORWARD=TRUE
    BACKWARD=FALSE
    var_in_mod<-1
    #fm0 <- paste(event, "~", paste0(var_in_mod, collapse = "+"),  sep = "")
    #mod0<-glm(fm0, data=data0, family=binomial)
  }
  if(direction%in%c("back_for","backward")){
    FORWARD=FALSE
    BACKWARD=TRUE
    var_in_mod<-varname
    #fm0 <- paste(event, "~", paste0(var_in_mod, collapse = "+"),  sep = "")
    #mod0<-glm(fm0, data=data0, family=binomial)
  }
  
  # Output matrix
  outmat<-data.frame(matrix(NA, ncol = 5, nrow=1))
  names(outmat)<-c("nstep","Var_in_model","VarIn","VarOut", "p")
  
  nstep<-1
  var_in_mod0<-var_in_mod
  while(TRUE){
    # Compare var_in_model before and after last step. If same, then break
    if(nstep!=1 & 
       length(setdiff(var_in_mod0, var_in_mod))==0 &
       length(setdiff(var_in_mod, var_in_mod0))==0){
      break
    }else{
      var_in_mod0<-var_in_mod
    }
    
    if(FORWARD==TRUE){
      var_list<-varname[varname%in%var_in_mod==0]
      var_p<-data.frame(matrix(NA,ncol=2,nrow=length(var_list)))
      names(var_p)<-c("VarName","p")
      
      # Fit modl w/ each of addition var and select the most sig
      for(i in 1:length(var_list)){
        var_i<-var_list[i]
        fm_i<-paste(event, "~", paste0(var_in_mod, collapse = "+"), "+",var_i,  sep = "")
        mod_i<-glm(fm_i, data=data0, family=binomial)
        tyoe3aov<-Anova(mod_i, test="Wald", type="III")
        var_p$VarName[i]<-rownames(tyoe3aov)[length(var_in_mod)+1]
        var_p$p[i]<-tyoe3aov$`Pr(>Chisq)`[length(var_in_mod)+1]
      } # end for loop
      p_in<-var_p$p[var_p$p==min(var_p$p)]
      var_in<-var_p$VarName[var_p$p==min(var_p$p)&var_p$VarName!="(Intercept)"]
      if(p_in<alpha_in){
        out_p<-ifelse(p_in<0.001, "<0.001",
                      as.character(round(p_in,3)))
        var_in_mod<-c(var_in_mod, var_in)
        outvec<-list(nstep, paste0(var_in_mod,collapse="+"), var_in, NA, out_p)
        outmat<-rbind(outmat, unlist(outvec))
        
        FORWARD=TRUE
        
        if(length(var_in_mod)==(ncol(data0)-1)){
          if(direction=="backward"){
            break
          }else if(direction%in%c("for_back","back_for")){
            # Change to backward
            FORWARD=TRUE
            BACKWARD=FALSE
          }}
      }else {
        if(direction=="forward"){
          break
        }else if(direction%in%c("for_back","back_for")){
            # Change to backward
            FORWARD=FALSE
            BACKWARD=TRUE
        }
      }
    } # end if(FORWARD=T)
    if(BACKWARD==TRUE){
      var_list<-var_in_mod
      var_p<-data.frame(matrix(NA,ncol=2,nrow=length(var_list)))
      names(var_p)<-c("VarName","p")
      
      fm_i <- paste(event, "~", paste0(var_list, collapse = "+"),  sep = "")
      mod_i<-glm(fm_i, data = data0, family=binomial)
      tyoe3aov<-Anova(mod_i, test="Wald", type="III")[-1,]
      if(min(tyoe3aov$`Pr(>Chisq)`==1)){ # if all se inflated but still converge, then break
        break
      }
      p_out<-tyoe3aov$`Pr(>Chisq)`[tyoe3aov$`Pr(>Chisq)`==max(tyoe3aov$`Pr(>Chisq)`)]
      var_out<-rownames(tyoe3aov)[tyoe3aov$`Pr(>Chisq)`==max(tyoe3aov$`Pr(>Chisq)`)]
      
      if(p_out>alpha_out){
        out_p<-ifelse(p_out<0.001, "<0.001",
                      as.character(round(p_out,3)))
        var_in_mod<-var_in_mod[-which(var_in_mod==var_out)]
        
        outvec<-list(nstep, paste0(var_in_mod,collapse="+"), NA, var_out, out_p)
        outmat<-rbind(outmat, unlist(outvec))
        
        if(length(var_in_mod)==0){
          if(direction=="backward"){
            break
          }else if(direction%in%c("for_back","back_for")){
            # Change to backward
            FORWARD=TRUE
            BACKWARD=FALSE
          }
        }
        
      }else {
        if(direction=="backward"){
          break
        }else if(direction%in%c("for_back","back_for")){
          # Change to backward
          FORWARD=TRUE
          BACKWARD=FALSE
        }
      }
      
    } # end if(BACKWARD=T)
    
    nstep<-nstep+1
    
  } # end while
  outmat<-outmat[-1,]
  return(outmat)

}
rcv<-function(datax, number, repeats, method, formulax, event, n_met, perform_metx, 
              nthreadx, ...){
  metrics_mat<-matrix(NA,nrow=repeats*number, ncol=n_met)
  model_rcv<-list(NA)
  
  for(r in 1:repeats){ # Repeat CV for xxx times
    grp<-dismo::kfold(c(1:nrow(datax)), k=number) # randomly split data into kfold
    fold_grp<-data.frame(cbind(c(1:nrow(datax)), grp))
    names(fold_grp)<-c("id","grp")
    for(f in 1:number){ # k-fold CV
      train_dat<-datax[-c(fold_grp$id[fold_grp$grp==f]),]
      test_dat<-datax[c(fold_grp$id[fold_grp$grp==f]),]
      
      if(method=="rfsrc"){
        m1<-rfsrc(formulax, data=train_dat, ...)
        p1<-predict(m1, newdata = test_dat, na.action = "na.impute")
        metrics_mat[(number*(r-1)+f),]<-perform_metx(y=p1$yvar, yhat=p1$class, phat=p1$predicted[,1])
        model_rcv[[number*(r-1)+f]]<-m1
      }
      if(method=="xgb"){
        # https://github.com/topepo/caret/blob/master/models/files/xgbTree.R
        out_vec<-ifelse(train_dat[,event]=="No",0,
                        ifelse(train_dat[,event]=="Yes",1,NA))
        cov_mat_train<-train_dat[,-which(colnames(train_dat)%in%c(event))]
        cov_mat_test<-test_dat[,-which(colnames(test_dat)%in%c(event))]
        
        # Build model
        if(is.factor(cov_mat_train)==1){
          ff<-~.
          mf<-model.frame(formula = ff, data = cov_mat_train, na.action = "na.pass")
          cov_mat_v1<-model.matrix(object = ff, data = mf)
          
          xgb_v1 <- xgboost(data = cov_mat_v1, label = out_vec, verbose = F, missing = NA, 
                            objective = "binary:logistic", na.action = "na.pass", nthreadx=nthreadx, ...)
          
          # Predict on test
          mf2 <- model.frame(formula = ff, data = cov_mat_test, na.action = "na.pass")
          cov_mat_v2 <- model.matrix(object = ff, data = mf2)
          test_dat$phat<-predict(xgb_v1, newdata = cov_mat_v2, missing = NA)
          test_dat$yhat<-factor(as.character(ifelse(test_dat$phat>=0.5,"Yes","No")), levels = c("No", "Yes"))
          metrics_mat[(number*(r-1)+f),]<-perform_metx(y=test_dat[,event], yhat=test_dat$yhat, phat=test_dat$phat)
          model_rcv[[number*(r-1)+f]]<-xgb_v1
        }else{
          cov_mat_v1 <- data.matrix(cov_mat_train)
          xgb_v1 <- xgboost(data = cov_mat_v1, label = out_vec, verbose = F, missing = NA, 
                            objective = "binary:logistic", na.action = "na.pass", nthreadx=nthreadx, ...)
          
          # Predict on test
          #mf2 <- model.frame(formula = ff, data = cov_mat_test, na.action = "na.pass")
          cov_mat_v2 <- data.matrix(cov_mat_test)
          test_dat$phat<-predict(xgb_v1, newdata = cov_mat_v2, missing = NA)
          test_dat$yhat<-factor(as.character(ifelse(test_dat$phat>=0.5,"Yes","No")), levels = c("No", "Yes"))
          metrics_mat[(number*(r-1)+f),]<-perform_metx(y=test_dat[,event], yhat=test_dat$yhat, phat=test_dat$phat)
          model_rcv[[number*(r-1)+f]]<-xgb_v1
        }
      }
      if(method=="BART"){
        cov_mat<-data.frame(train_dat[,-which(names(train_dat)==event)])
        names(cov_mat)<-names(train_dat)[-which(names(train_dat)==event)]
        train_dat[,event]<-factor(as.character(train_dat[,event]), levels = c("Yes","No"))
        bm1<-bartMachine(y = train_dat[,event], X=cov_mat, 
                         mem_cache_for_speed=FALSE,
                         use_missing_data = TRUE, use_missing_data_dummies_as_covars = TRUE) # 1min
        
        cov_mat_test<-data.frame(test_dat[,-which(names(test_dat)==event)])
        names(cov_mat_test)<-names(test_dat)[-which(names(test_dat)==event)]
        
        p1<-predict(bm1, new_data=cov_mat_test)
        c1<-predict(bm1, new_data=cov_mat_test, type="class")
        
        metrics_mat[(number*(r-1)+f),]<-perform_metx(y=test_dat[,event], yhat=c1, 
                                                     phat=p1)
        model_rcv[[number*(r-1)+f]]<-bm1
        
        
      }
      if(method=="glm"){
        m1<-glm(formulax, data=train_dat, family = binomial)
        test_dat$phat<-predict(m1, newdata=test_dat, type="response")
        test_dat$yhat<-factor(as.character(ifelse(test_dat$phat>=0.5,"Yes","No")), levels = c("No", "Yes"))
        metrics_mat[(number*(r-1)+f),]<-perform_metx(y=test_dat[,event], yhat=test_dat$yhat, phat=test_dat$phat)
        model_rcv[[number*(r-1)+f]]<-m1
      }
      
      cat(paste("Rep ",r," Fold ",f,"\n",sep=""))
    }
  }
  colnames(metrics_mat)<-colnames(perform_metx(y=factor(c("No","Yes")), yhat=factor(c("No","Yes")), phat=c(0.3,0.2)))
  rownames(metrics_mat)<-paste(paste("Rep",rep(c(1:repeats),each=number),sep="_"), paste("Fold",c(1:number),sep="_"), sep="_")
  names(model_rcv)<-paste(paste("Rep",rep(c(1:repeats),each=number),sep="_"), paste("Fold",c(1:number),sep="_"), sep="_")
  
  return(list("model_rcv"=model_rcv, "metrics_mat"=metrics_mat))
}
var_select_backward<-function(datax, event, method, number=10, repeats=1, var_drop_rate=0.1,
                              mtry_form=ceiling(sqrt(nvar)), nroundx=200, c.sd=1, rcv_ce=FALSE, nthreadx=1,
                              cf_threshold=NULL, ntreex1=5000, ntreex2=2000, minc_x=NULL){
  datax1<-datax
  datax1$event<-datax1[,c(event)]
  datax1<-datax1[,-which(names(datax1)==event)]
  
  metrics_mat<-data.frame(matrix(NA,nrow=ncol(datax1)-1, ncol=1))
  names(metrics_mat)<-"met"
  model_list<-list(NA)
  #imp_list<-list(NA)
  var_out_list<-list(NA)
  
  if(method=="xgb"){
    out_vec<-ifelse(datax1$event=="No",0,
                    ifelse(datax1$event=="Yes",1,NA))
    
    # Build a full model to get var importance
    cov_mat0<-datax1[,-which(colnames(datax1)%in%c("event"))]
    cov_mat_v0<-data.matrix(cov_mat0)
    xgb_v0 <- xgboost(data = cov_mat_v0, label = out_vec, verbose = F, missing = NA, 
                      objective = "binary:logistic", na.action = "na.pass", 
                      nrounds = nroundx, nthread=nthreadx)
    
    imp_list0<-xgb.importance(model=xgb_v0)
    if(nrow(imp_list0)<ncol(cov_mat0) & round(sum(imp_list0$Frequency))==1){
      var_not_in_imp<-names(cov_mat0)[(names(cov_mat0)) %in% imp_list0$Feature==0]
      mock_mat<-data.frame(matrix(NA, nrow=length(var_not_in_imp), ncol=4))
      names(mock_mat)<-colnames(imp_list0)
      mock_mat[,1]<-var_not_in_imp
      mock_mat[,c(2:4)]<-rep(0,length(var_not_in_imp)*3)
      imp_list0<-rbind(imp_list0, mock_mat)
    }
    imp_list<-imp_list0
    
    imp_list_i<-list(imp_list)
    i<-1
    nvar<-ncol(datax1)-1
    
    # Random split data into training and test for xgb (1:1)
    sam_id<-sample(c(1:nrow(datax1)), round(nrow(datax1)/2))
    dat_train<-datax1[sam_id,]
    dat_test<-datax1[-sam_id,]
    out_vec_train<-ifelse(dat_train$event=="No",0,
                          ifelse(dat_train$event=="Yes",1,NA))
    
    # Var selection
    while(nvar>=1){ 
      cov_mat<-datax1[,-which(colnames(datax1)%in%c("event"))]
      cov_mat_train<-dat_train[,-which(colnames(dat_train)%in%c("event"))]
      cov_mat_test<-dat_test[,-which(colnames(dat_test)%in%c("event"))]
      
      metrics_mat$var_in_mod[i]<-paste0(colnames(datax1)[colnames(datax1)%in%c("event")==0],collapse = "_")
      metrics_mat$numcov[i]<-nvar
      
      if(is.factor(cov_mat)==1){ # If only a factor var left, model can't run
        
        # Calculate classification error from 10-fold CV
        if(rcv_ce==TRUE){
          xgb_cv<-rcv(datax=datax1, number=number, repeats=repeats, method="xgb", formulax=event~., event = "event",
                      perform_metx=class_error, n_met=1, nrounds=nroundx, nthreadx=nthreadx)
          metric<-data.frame(xgb_cv[[2]])
          metrics_mat$met[i]<-apply(metric,2,mean)
        }else{
          if(is.null(cov_mat_train)==1){
            cov_mat_train<-data.frame(cov_mat_train)
            cov_mat_test<-data.frame(cov_mat_test)
            names(cov_mat_train)<-names(cov_mat_test)<-names(datax1)[names(datax1)!="event"]
          }
          
          # Fit using train
          ff<-~.
          mf_train<-model.frame(formula = ff, data = data.frame(cov_mat_train), na.action = "na.pass")
          names(mf_train)<-names(datax1)[names(datax1)!="event"]
          cov_mat_v1_train<-model.matrix(object = ff, data = mf_train)
          
          mf_test<-model.frame(formula = ff, data = data.frame(cov_mat_test), na.action = "na.pass")
          names(mf_test)<-names(datax1)[names(datax1)!="event"]
          cov_mat_v1_test<-model.matrix(object = ff, data = mf_test)
          
          xgb_v1 <- xgboost(data = cov_mat_v1_train, label = out_vec_train, verbose = F, missing = NA, 
                            objective = "binary:logistic", na.action = "na.pass", nrounds=nroundx, nthreadx=nthreadx)
          #model_list[[i]]<-xgb_v1
          
          # Predict on test
          phat<-predict(xgb_v1, newdata = cov_mat_v1_test, missing = NA)
          yhat<-factor(as.character(ifelse(phat>=0.5,"Yes","No")), levels = c("No", "Yes"))
          
          metrics_mat$met[i]<-class_error(y=dat_test[,"event"], yhat=yhat, phat=phat)
          
        } # end else
        
        # Var imp
        var_out<-"0"
        var_out_list[[i]]<-var_out
        
        # Exclude var from data
        datax1<-datax1[,-which(names(datax1)%in%var_out)]
      }else{  # else for (if is.factor(cov_mat)==1)
        if(rcv_ce==TRUE){
          # Calculate classification error  from 10-fold CV
          xgb_cv<-rcv(datax=datax1, number=number, repeats=repeats, method="xgb", formulax=event~., event = "event",
                      perform_metx=class_error, n_met=1, nrounds=nroundx, nthreadx=nthreadx)
          metric<-data.frame(xgb_cv[[2]])
          metrics_mat$met[i]<-apply(metric,2,mean)
        }else{
          # Fit using train
          if(is.null(ncol(cov_mat_train))==1){
            cov_mat_train<-data.frame(cov_mat_train)
            cov_mat_test<-data.frame(cov_mat_test)
            names(cov_mat_train)<-names(cov_mat_test)<-names(datax1)[names(datax1)!="event"]
          }
          
          ff<-~.
          mf_train<-model.frame(formula = ff, data = data.frame(cov_mat_train), na.action = "na.pass")
          cov_mat_v1_train<-model.matrix(object = ff, data = mf_train)
          
          mf_test<-model.frame(formula = ff, data =  data.frame(cov_mat_test), na.action = "na.pass")
          cov_mat_v1_test<-model.matrix(object = ff, data = mf_test)
          
          xgb_v1 <- xgboost(data = cov_mat_v1_train, label = out_vec_train, verbose = F, missing = NA, 
                            objective = "binary:logistic", na.action = "na.pass", nrounds=nroundx, nthreadx=nthreadx)
          
          # Predict on test
          phat<-predict(xgb_v1, newdata = cov_mat_v1_test, missing = NA)
          yhat<-factor(as.character(ifelse(phat>=0.5,"Yes","No")), levels = c("No", "Yes"))
          
          metrics_mat$met[i]<-class_error(y=dat_test[,"event"], yhat=yhat, phat=phat)
        } # end else
        #model_list[[i]]<-xgb_v1
        
        
        # Var imp
        if(is.null(ncol(cov_mat))==0){
          num_var_exclude<-round(var_drop_rate*ncol(cov_mat))
          if(num_var_exclude==0){num_var_exclude=1}
          imp_list<-imp_list[order(imp_list$Gain, decreasing=FALSE),]
          var_out<-imp_list$Feature[c(1:num_var_exclude)]
          var_out_list[[i]]<-var_out
        }else{
          num_var_exclude<-1
          imp_list<-imp_list[order(imp_list$Gain, decreasing=FALSE),]
          var_out<-imp_list$Feature[c(1:num_var_exclude)]
          var_out_list[[i]]<-var_out
        }
        
        
        # Exclude var from data
        datax1<-datax1[,-which(names(datax1)%in%var_out)]
        dat_train<-dat_train[,-which(names(dat_train)%in%var_out)]
        dat_test<-dat_test[,-which(names(dat_test)%in%var_out)]
        
        imp_list_i[[i+1]]<-imp_list<-imp_list[imp_list$Feature%in%var_out==0,]
      }
      cat(paste(" i ",i,"\n",sep=""))
      i<-i+1
      nvar<-nvar-num_var_exclude
    }
    
    met_dat<-na.omit(metrics_mat)
    met_dat$met_sd<-sqrt(met_dat$met*(1-met_dat$met)/nrow(datax))
    min_met<-min(met_dat$met)
    min_met_ci<-min_met+c.sd*met_dat$met_sd[which.min(met_dat$met)]
    best_var_list_row<-which(met_dat$met <= min_met_ci)[which.min(met_dat$numcov[which(met_dat$met <= 
                                                                                         min_met_ci)])]
    best_var_list<-met_dat$var_in_mod[best_var_list_row]
    var_in_model<-strsplit(best_var_list,"_")
    
  }
  if(method=="cforest"){
    mtryx<-mtry_form(nrow(datax1))
    cf0<-cforest(event~., data=datax1, control=ctree_control(MIA = FALSE), 
                 ntree=ntreex1, mtry = mtryx)
    #cf0<-cforest(event~., data=datax1, control=ctree_control(MIA = TRUE))
    #imp_list0<-varimp(cf0, conditional = TRUE, threshold = cf_threshold)
    imp_list0<-varimp(cf0, conditional = TRUE)
    imp_list0<-imp_list0[names(imp_list0)!=""]
    if(length(imp_list0)<(ncol(datax1)-1)){
      var_not_in_imp<-names(datax1)[(names(datax1)) %in% c("event",names(imp_list0))==0]
      if(min(imp_list0)<=0){
        mock_mat<-rep((min(imp_list0)-1), length(var_not_in_imp))
      }else{
        mock_mat<-rep(0, length(var_not_in_imp))
      }
      names(mock_mat)<-var_not_in_imp
      imp_list0<-c(imp_list0, mock_mat)
    }
    
    imp_list<-imp_list0
    imp_list_i<-list(imp_list)
    i<-1
    nvar<-ncol(datax1)-1
    while(nvar>=1){ 
      metrics_mat$var_in_mod[i]<-paste0(colnames(datax1)[colnames(datax1)%in%c("event")==0],collapse = "_")
      metrics_mat$numcov[i]<-nvar
      #mtryx<-mtry_form(nrow(datax1))
      
      if(i==1){
        cf1<-cf0
      }else{
        cf1<-cforest(event~., data=datax1, control=ctree_control(MIA = FALSE), 
                     ntree=ntreex2, mtry = mtryx)
        #cf1<-cforest(event~., data=datax1, control=ctree_control(MIA = FALSE))
      }
      oobpredict<-predict(cf1, OOB = TRUE)
      oob_ce<-class_error(y=datax1$event, yhat=oobpredict)
      metrics_mat$met[i]<-oob_ce
      model_list[[i]]<-cf1
      
      # Var imp
      if(is.null(ncol(datax1))==0 | ncol(datax1)!=1){
        num_var_exclude<-round(var_drop_rate*(ncol(datax1)-1))
        if(num_var_exclude==0){num_var_exclude=1}
        if(sum(imp_list<=0)<=num_var_exclude){
          imp_list<-imp_list[order(imp_list, decreasing=FALSE)]
          var_out<-names(imp_list)[c(1:num_var_exclude)]
          var_out_list[[i]]<-var_out
        }else{
          imp_list<-imp_list[order(imp_list, decreasing=FALSE)]
          var_out<-names(imp_list)[imp_list<=0]
          var_out_list[[i]]<-var_out
        }
        
        
      }else{
        num_var_exclude<-1
        imp_list<-imp_list[order(imp_list, decreasing=FALSE)]
        var_out<-names(imp_list)[c(1:num_var_exclude)]
        var_out_list[[i]]<-var_out
      }
      
      # Exclude var from data
      num_var_exclude<-length(var_out)
      datax1<-datax1[,-which(names(datax1)%in%var_out)]
      imp_list_i[[i+1]]<-imp_list<-imp_list[names(imp_list)%in%var_out==0]
      
      cat(paste(" i ",i,"\n",sep=""))
      i<-i+1
      nvar<-nvar-num_var_exclude
    }
    met_dat<-na.omit(metrics_mat)
    met_dat$met_sd<-sqrt(met_dat$met*(1-met_dat$met)/nrow(datax))
    min_met<-min(met_dat$met)
    min_met_ci<-min_met+c.sd*met_dat$met_sd[which.min(met_dat$met)]
    best_var_list_row<-which(met_dat$met <= min_met_ci)[which.min(met_dat$numcov[which(met_dat$met <= 
                                                                                         min_met_ci)])]
    best_var_list<-met_dat$var_in_mod[best_var_list_row]
    var_in_model<-strsplit(best_var_list,"_")
    
    
  }
  return(list("metrics_mat"=met_dat, "varimp"=imp_list0, 
              "var_in_model"=var_in_model))
  #"model_list"=model_list
}
boot_finalvar<-function(cutpt_boot_varselect, varlist, n_boot){
  final_var<-list(NA)
  met_max<-list(NA)
  
  for(j in 1:length(cutpt_boot_varselect)){
    boot_cut_j<-cutpt_boot_varselect[j]
    varlist_i<-table(unlist(varlist))
    # Identify final var by cut_j
    final_var[[j]]<-names(varlist_i[varlist_i>=boot_cut_j*n_boot])
  }
  return(final_var)
}


bisectK <- function(tol, coverage, permute_mat, x_left, x_right, countLimit, perm_mean, perm_se){
  count = 0
  guess = mean(c(x_left, x_right))
  while ((x_right - x_left) / 2 >= tol & count < countLimit){
    empirical_coverage = mean(sapply(1 : nrow(permute_mat), 
                                     function(s){all(permute_mat[s,] - perm_mean <= guess * perm_se)}))
    if (empirical_coverage - coverage == 0){
      break
    } else if (empirical_coverage - coverage < 0){
      x_left = guess
    } else {
      x_right = guess
    }
    guess = mean(c(x_left, x_right))
    count = count + 1
  }
  guess
}

# Metrics
prec_recall_f1<-function(gold_pos,gold_neg,var_selected){
  tp<-intersect(gold_pos,var_selected); tp_n<-length(tp)
  fp<-intersect(gold_neg,var_selected); fp_n<-length(fp)
  fn<-setdiff(gold_pos,var_selected); fn_n<-length(fn)
  precision<-tp_n/(tp_n+fp_n)
  recall<-tp_n/(tp_n+fn_n)
  f1<-2*precision*recall/(precision+recall)
  return(list("tp"=tp,"fp"=fp,"fn"=fn,"precision"=precision,"recall"=recall,"f1"=f1))
}
