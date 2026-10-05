# Run from the SEM directory after obtaining authorized access to the analysis dataset.
# Rscript reproduce_sem.R /path/ChpSTREM2AD_analysis_dataset.csv /path/output [bootstrap]
# 从 SEM 文件夹运行，models 中的文件名用于识别各条路径。
# モデル内の矢印は横断データの関連を表すもので、時間的順序を確定するものではない。
if(.Platform$OS.type=="windows")invisible(try(Sys.setlocale("LC_CTYPE","English_United States.utf8"),silent=TRUE))
args<-commandArgs(trailingOnly=TRUE)
if(length(args)<2)stop("Provide the harmonized dataset CSV and an output directory.")
if(nzchar(Sys.getenv("SEM_R_LIBRARY"))).libPaths(c(Sys.getenv("SEM_R_LIBRARY"),.libPaths()))
suppressPackageStartupMessages(library(lavaan))
source_dir<-getwd();deep<-normalizePath(args[2],mustWork=FALSE)
dir.create(deep,recursive=TRUE,showWarnings=FALSE)
d<-read.csv(args[1],check.names=FALSE)
required<-c("MMSE","MOCA","ADAS13","S_PTGENDER","S_AGE","L_PTEDUCAT","APOE401","L_ChPICV","L_MSD_STREM2CORRECTED","L_TAU","L_PTAU","L_S_ABETA","S_DX_label")
stopifnot(all(required %in% names(d)),nrow(d)==735)
# L_ 列已经过对数转换，这里不再重复取 log。
lookup<-c(ChP="L_ChPICV",sTREM2="L_MSD_STREM2CORRECTED",TAU="L_TAU",PTAU="L_PTAU",ABETA="L_S_ABETA")
files<-list.files(file.path(source_dir,"models"),pattern="\\.lavaan$",full.names=TRUE)
estimates<-list();diagnostics<-list();jobs<-list()
for(file in files) {
  id<-sub("\\.lavaan$","",basename(file));parts<-strsplit(id,"__",fixed=TRUE)[[1]];kind<-parts[2];model<-parts[3];gr<-parts[4]
  v<-strsplit(model,"_to_",fixed=TRUE)[[1]]
  if(kind=="single")input<-data.frame(x=d[[lookup[[v[1]]]]],mediator=d[[lookup[[v[2]]]]])
  if(kind=="parallel")input<-data.frame(x=d$L_ChPICV,mediator_1=d$L_MSD_STREM2CORRECTED,mediator_2=d[[if(grepl("PTAU",model))"L_PTAU" else "L_TAU"]],mediator_3=d$L_S_ABETA)
  if(kind=="serial")input<-data.frame(x=d[[lookup[[v[1]]]]],mediator_1=d$L_MSD_STREM2CORRECTED,mediator_2=d$L_ChPICV)
# ADAS13 は高いほど成績が悪いので符号を反転し、他の2指標と向きをそろえる。
  input$indicator_1<-d$MMSE;input$indicator_2<-d$MOCA;input$indicator_3<--d$ADAS13
  input$S_PTGENDER<-d$S_PTGENDER;input$S_AGE<-d$S_AGE;input$PTEDUCAT<-d$L_PTEDUCAT;input$APOE401<-d$APOE401
# 标准化在全样本中完成，随后才加入分组标签；不是各诊断组分别标准化。
  input<-as.data.frame(lapply(input,function(x)as.numeric(scale(x))))
  model_text<-paste(readLines(file,warn=FALSE),collapse="\n")
# 点推定はMLR、欠測はFIMLで扱う。bootstrap内の推定設定はbootstrap_fit.Rを参照。
  a<-list(model=model_text,data=input,estimator="MLR",missing="fiml",fixed.x=TRUE,std.lv=TRUE,meanstructure=TRUE)
# 多群モデルでは因子負荷量を群間で等しくする。構造パスまで一律に固定しない。
  if(gr=="Multigroup") {input$.group<-d$S_DX_label;a$data<-input;a$group<-".group";a$group.label<-c("CN","MCI","AD");a$group.equal<-"loadings"}
  f<-do.call(sem,a);p<-parameterEstimates(f,standardized=TRUE);pp<-p;pp$id<-id;estimates[[id]]<-pp
  fm<-fitMeasures(f,c("df","cfi.robust","rmsea.robust","srmr"));diagnostics[[id]]<-data.frame(id=id,n=sum(lavInspect(f,"nobs")),converged=lavInspect(f,"converged"),t(fm),check.names=FALSE)
  jobs[[id]]<-list(id=id,name=model,kind=kind,definition=parts[1],group=gr,group_equal=if(gr=="Multigroup")"loadings" else NULL,data=input,model=model_text,point=p,start=coef(f))
}
# parameterEstimates has a group column only for multigroup models.
cols<-unique(unlist(lapply(estimates,names)))
estimates<-lapply(estimates,function(p){for(n in setdiff(cols,names(p)))p[[n]]<-if(n=="group")1 else NA;p[,cols]})
write.csv(do.call(rbind,estimates),file.path(deep,"point_estimates.csv"),row.names=FALSE)
write.csv(do.call(rbind,diagnostics),file.path(deep,"fit_indices.csv"),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(deep,"sessionInfo.txt"))
# bootstrap 耗时较长，只有明确传入第3个参数 bootstrap 才运行。
if(length(args)>=3 && args[3]=="bootstrap") {
  source(file.path(source_dir,"bootstrap_fit.R"),encoding="UTF-8")
  dir.create(file.path(deep,"bootstrap"),showWarnings=FALSE)
  status<-lapply(jobs,fit_job)
  write.csv(do.call(rbind,status),file.path(deep,"bootstrap_status.csv"),row.names=FALSE)
}
cat("SEM reproduction complete\n")
