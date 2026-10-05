# 每次重抽样都重新拟合模型，并保留收敛和异常方差等诊断。
# 既存CSVがある場合は途中から再開する。同じ出力先を使うときは入力とモデルも同じであることを確認。
fit_job<-function(j) {
  suppressPackageStartupMessages(library(lavaan))
  B<-1000L
  target<-file.path(deep,"bootstrap",paste0(j$id,".csv"))
  if(file.exists(target)) {
    prev<-read.csv(target,check.names=FALSE)
    if(nrow(prev)==B) return(data.frame(id=j$id,n=B,valid=sum(prev$valid),negative=sum(prev$negative_variances>0,na.rm=TRUE),nonconverged=sum(!prev$converged)))
  } else prev<-NULL
  set.seed(20260926L + match(j$group,c("Overall","CN","MCI","AD","Multigroup"))*1000L)
  n<-nrow(j$data)
# 多群モデルは診断群内で復元抽出し、各群の人数を維持する。
  indices<-if(j$group=="Multigroup")replicate(B,unlist(lapply(split(seq_len(n),j$data$.group),function(ii)sample(ii,length(ii),replace=TRUE)),use.names=FALSE),simplify=FALSE) else replicate(B,sample.int(n,n,replace=TRUE),simplify=FALSE)
  extract<-function(pe) {
    def<-pe[pe$op==":=",]
    ans<-c(setNames(def$est,paste0(def$lhs,"__est")),setNames(def$std.all,paste0(def$lhs,"__std")))
    pth<-pe[pe$op=="~" & pe$label!="",]
    ans<-c(ans,setNames(pth$est,paste0(pth$label,"__est")),setNames(pth$std.all,paste0(pth$label,"__std")))
    lam<-pe[pe$op=="=~",]
    ln<-if("group" %in% names(lam) && length(unique(lam$group))>1)paste0(lam$rhs,"_g",lam$group,"__loading") else paste0(lam$rhs,"__loading")
    c(ans,setNames(lam$std.all,ln))
  }
  fields<-names(extract(j$point));rows<-vector("list",B)
  if(!is.null(prev)) for(k in seq_len(nrow(prev))) rows[[k]]<-prev[k,,drop=FALSE]
  first<-if(is.null(prev))1L else nrow(prev)+1L
  for(b in seq.int(first,B)) {
    warns<-character();err<-""
    fitit<-function(initial) {
      args<-list(model=j$model,data=j$data[indices[[b]],,drop=FALSE],estimator="ML",missing="fiml",fixed.x=TRUE,std.lv=TRUE,meanstructure=TRUE,se="none",test="none")
      if(j$group=="Multigroup") {args$group<-".group";args$group.label<-c("CN","MCI","AD");if(!is.null(j$group_equal))args$group.equal<-j$group_equal}
      if(initial)args$start<-j$start
      tryCatch(withCallingHandlers(do.call(lavaan::sem,args),warning=function(w){warns<<-c(warns,conditionMessage(w));invokeRestart("muffleWarning")}),error=function(e){err<<-conditionMessage(e);NULL})
    }
# 先尝试原模型的起始值；未收敛时再用默认起始值重试一次。
    f<-fitit(TRUE);retry<-FALSE
    if(is.null(f) || !isTRUE(lavInspect(f,"converged"))) {retry<-TRUE;f<-fitit(FALSE)}
    vals<-setNames(rep(NA_real_,length(fields)),fields)
    conv<-FALSE;negative<-NA_integer_;grad<-NA_real_;rank<-NA_integer_;np<-NA_integer_;orientation<-FALSE
    if(!is.null(f)) {
      conv<-isTRUE(lavInspect(f,"converged"))
      pe<-suppressWarnings(parameterEstimates(f,standardized=TRUE)); vv<-pe[pe$op=="~~" & pe$lhs==pe$rhs,]
      negative<-sum(vv$est < -1e-8,na.rm=TRUE)
      score<-lavInspect(f,"gradient")
# 等値制約がある場合、勾配は制約を満たす方向に射影して確認する。
      cj<-f@Model@con.jac
      if(is.matrix(cj) && nrow(cj)>0 && ncol(cj)==length(score)) {
        qq<-qr(t(cj));q<-qr.Q(qq,complete=TRUE)
        score<-as.numeric(crossprod(q[,seq.int(qq$rank+1L,ncol(q)),drop=FALSE],score))
      }
      grad<-max(abs(score))
# 因子负荷保持正向，避免潜变量符号翻转后直接混合汇总。
      orientation<-all(pe$est[pe$op=="=~"]>0)
      delta<-do.call(rbind,lavTech(f,"delta"))
      if(is.matrix(cj) && nrow(cj)>0 && ncol(cj)==ncol(delta))delta<-rbind(delta,cj)
      cs<-sqrt(colSums(delta^2));cs[cs==0]<-1
      ds<-svd(sweep(delta,2,cs,"/"),nu=0,nv=0)$d;rank<-sum(ds>max(ds)*1e-8);np<-ncol(delta)
      v<-extract(pe);vals[names(v)]<-v
    }
# 收敛不等于可用；还要检查方差、梯度、局部秩和参数是否有限。
    valid<-conv && isTRUE(negative==0) && isTRUE(grad<1e-4) && isTRUE(rank==np) && orientation && all(is.finite(vals))
    reason<-if(!conv)"nonconvergence" else if(negative>0)"negative_variance" else if(grad>=1e-4)"gradient" else if(rank<np)"rank" else if(!orientation)"orientation" else if(!all(is.finite(vals)))"nonfinite" else "ok"
    rows[[b]]<-data.frame(iteration=b,valid=valid,converged=conv,negative_variances=negative,max_gradient=grad,rank=rank,npar=np,retried=retry,reason=reason,warning_count=length(unique(warns)),error=err,as.list(vals),check.names=FALSE)
# 50回ごとに保存。中断時にも計算済みの反復を残す。
    if(b%%50==0 || b==B) {
      write.csv(do.call(rbind,rows[seq_len(b)]),target,row.names=FALSE,na="")
      cat(j$id,b,"/",B,"\n");flush.console()
    }
  }
  rr<-do.call(rbind,rows)
  data.frame(id=j$id,n=B,valid=sum(rr$valid),negative=sum(rr$negative_variances>0,na.rm=TRUE),nonconverged=sum(!rr$converged))
}
