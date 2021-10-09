
nameToPeopleLowHigh=function(n) {
  if(n=="Arches" || n=="Dock") {
    return(1)
  } else {
    return(0)
  }
}


addRobustSEToLMModel=function(lmm) {
lmm$vcovHC_ <- vcovHC(lmm)
lmm$robust.se <-  sqrt(diag(lmm$vcovHC_))
lmm$vcovCONSTANT_=vcovHC(lmm, "const")
lmm$seCONSTANT_=sqrt(diag(lmm$vcovCONSTANT_))
return(lmm)
}


createRegressionData=function(dtable) {
  names=c("Arches","Beach","Dock","Rainforest")
  #dtable=sortLvls.fnc(dtable$Name, c("Beach","Arches","Dock","Rainforest"))
  reg_data_people_scene=c()
  reg_data_people_indic=c()
  reg_data_people_low=c()
  reg_data_people_high=c()
  reg_data_people_click=c()
  for(n_i in 1:length(names)) {
    for(p_i in 0:1) {
      name_row=dtable[dtable$Name==names[n_i],]
      the_row=name_row[name_row$People==p_i,]
      #print("The row is ")
      #print(the_row)
      impressions=the_row$Impressions
      clicks=the_row$Results
      #hl_stat=
      print(paste(names[n_i],";",p_i," imp ",impressions," ; ",clicks," ; ",nameToPeopleLowHigh(names[n_i])))
      #print("new")
      for(impression in 1:impressions) {
        if(impression<=clicks) {
          reg_data_people_click=c(reg_data_people_click,1)
        } else {
          reg_data_people_click=c(reg_data_people_click,0)
        }
        reg_data_people_scene=c(reg_data_people_scene,names[n_i])
        reg_data_people_indic=c(reg_data_people_indic,p_i)
        reg_data_people_high=c(reg_data_people_high,nameToPeopleLowHigh(names[n_i]))
      }
    }
  }
  reg_df=data.frame(Scene=reg_data_people_scene,People=reg_data_people_indic,PeopleHigh=reg_data_people_high,Click=reg_data_people_click)
  #return(reg_df)
  write.csv(reg_df,"reg_data_df.csv", row.names = FALSE)
  
}


loadSavedRegDF=function() {
  reg_df=read.csv("reg_data_df.csv")
  #reg_df$Name=as.factor(reg_df$Name)
  #reg_df$Name=sortLvls.fnc(reg_df$Name, c("Beach","Arches","Dock","Rainforest"))
  return(reg_df)
}

getNamesFromDTable=function(dtable) {
  the_names=dtable$Name
  return(unique(the_names))
}

readAdCSV=function() {
   return(read.csv("adata.csv"))
}

getSubDFByName=function(n,bigdf) {
  #print(bigdf$Name==n)
  subdf=bigdf[bigdf$Name==n,]
  #print(subdf)
  return (subdf)
}

getNormSDFromData=function(imps,clicks) {
  p_a=clicks/imps
  samp_a=imps
  sqrt_num=(p_a)*(1.0-p_a)
  #print(paste("p_a",p_a))
  #print(paste("sqrt_num",sqrt_num))
  stddev_a=sqrt( sqrt_num / samp_a )
  #print(paste("stddev_a",stddev_a))
  return (stddev_a)
}

getImpsAndClicksFromRow=function(dfrow,impname) {
  clicks=dfrow$Results
  imps=dfrow[1,impname]
  return (c(imps,clicks))
}


calculateTwoSidedPValue=function(p_difference,diff_sd) {
  #p_difference is the difference between the two proportions
  #diff_sd is the standard deviation of the normal distribution
  #used to compute z-crit values.  
  #The mean of the distribution is zero.
    
  #pnorm gives area under the curve from negative infinity
  #up to the given point.  So to compute a two-sided p-value,
  #plug in the negative of the absolute value of p_difference, 
  #then, because it is two-sided, multiply it by two
  #print(paste("in2sp p_diff is ",p_difference,", mean=0 , diff_sd is ",diff_sd))
  pnorm_from_diff=pnorm(-abs(p_difference),mean=0,sd=abs(diff_sd))
  
  area_under_curve_to_the_left=pnorm_from_diff
  p_value=area_under_curve_to_the_left*2
  return(p_value)
}


generateTestText=function(np_click,np_imp,p_click,p_imp,name,first,alpha_thresh) {
  gen_text="<br><br>"

  np_rate=np_click/np_imp
  np_stddev=getNormSDFromData(np_imp,np_click)
  p_stddev=getNormSDFromData(p_imp,p_click)
  p_rate=p_click/p_imp
  gen_text=paste(gen_text,"For the ",name," set of Ads we have ",specify_decimal(np_click,0)," clicks and ",specify_decimal(np_imp,0)," ")
  gen_test=paste(gen_text," for no people for a click rate of ", specify_decimal(np_rate,4),".  ")
  gen_text=paste(gen_text,"For the ",name," set of Ads we have ",specify_decimal(p_click,0)," clicks and ",specify_decimal(p_imp,0)," ")
  gen_test=paste(gen_text," for people for a click rate of ", specify_decimal(p_rate,4),".  ")  
  if(first) {
    gen_text=paste(gen_text,"For proportions, the sampling distribution may be approximated with a normal distribution ",sep="")
    gen_text=paste(gen_text,"whose mean is the sampled proportion, but whose standard deviation $\\sigma$ is equal to ",sep="")
    gen_text=paste(gen_text,"","$\\sqrt{ \\frac{p(1-p)}{n}}$.  ",sep="")
  }
  gen_text=paste(gen_text,"  We have ",specify_decimal(np_stddev,4)," and ",specify_decimal(p_stddev,4)," as standard deviations for the sampling distributions for ",sep="")
  gen_text=paste(gen_text," no people and people respectively.  ",sep="")
                   
  #carry out prop.test
  successes=c(np_click,p_click)
  trials=c(np_imp,p_imp)
  prop_test_result=prop.test(successes,trials)
  prop_test_pval=prop_test_result$p.value
  if(first) {
    gen_text=paste(gen_text,"To compare the two proportions we may set $H_{0}:\\hat{p_{p}}=\\hat{p_{np}}$ .  ",sep="")
    gen_text=paste(gen_text,"We use the prop.test function to compute p-values in our testing.  ",sep="")
  }  
  if(prop_test_pval<alpha_thresh) {
    gen_text=paste(gen_text,"The p-value ",specify_decimal(prop_test_pval,4)," indicates that we may <em>reject</em> ",sep="")
  } else {
    gen_text=paste(gen_text,"The p-value ",specify_decimal(prop_test_pval,4)," indicates that we <em>fail</em> to reject ",sep="")
  }
  gen_text=paste(gen_text,"the null hypothesis that the two proportions are equal",sep="")
  gen_text=paste(gen_text,"at the ",sep="")
  if(first) {
    gen_text=paste(gen_text,"$\\alpha$=",sep="")
  }
  gen_text=paste(gen_text,specify_decimal(alpha_thresh,3)," confidence level.",sep="")
  return(c(gen_text,prop_test_pval))
}


generateTestSection=function(addata,name,isFirst) {
  colors=c("black","green")
  subdf=getSubDFByName(name,addata)
  people_row=subdf[subdf$People==1,]
  no_people_row=subdf[subdf$People==0,]
  people_imps_clicks=getImpsAndClicksFromRow(people_row,"Impressions")
  no_people_imps_clicks=getImpsAndClicksFromRow(no_people_row,"Impressions")
  people_sd=getNormSDFromData(people_imps_clicks[1],people_imps_clicks[2])
  no_people_sd=getNormSDFromData(no_people_imps_clicks[1],no_people_imps_clicks[2])
  people_rate=people_row$Click.Rate.Impressions
  no_people_rate=no_people_row$Click.Rate.Impressions
  x_data=seq(0,0.065,length=1000)
  y_vals_no_people=dnorm(x_data,mean=no_people_rate,sd=no_people_sd)
  y_vals_people=dnorm(x_data,mean=people_rate,sd=people_sd)
  plot_y_max=max(max(y_vals_no_people),max(y_vals_people))
  alpha=0.05
  dataTest=generateTestText(no_people_imps_clicks[2],no_people_imps_clicks[1],people_imps_clicks[2],people_imps_clicks[1],name,isFirst,alpha)
  ttp=dataTest[1]
  p_val=dataTest[2]
  plot(c(0), c(10), type="l", lty=2,ylim=c(0,plot_y_max*1.1),xlim=c(0,max(x_data)), xlab="Click Rate",ylab=paste("Density",name), main=paste("Comparison of",name,"Sampling Distributions"))
  lines(x_data,y_vals_no_people,col=colors[1])
  lines(x_data,y_vals_people,col=colors[2])
  leg_names=c(paste("Ctrl ; ",specify_decimal(no_people_rate,3)),paste("Trtmt ; ",specify_decimal(people_rate,3)))
  legend("topright", inset=.05, title="Distributions",leg_names, col=colors,lwd=1, lty=c(1, 1))
  cat(ttp)  
  return(as.double(p_val))
}


#https://stackoverflow.com/questions/3443687/formatting-decimal-places-in-r
specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))
