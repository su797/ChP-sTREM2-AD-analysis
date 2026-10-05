# 这里只读取汇总估计，不需要受试者级数据。
# 第1引数はaggregate_resultsを含むフォルダ。出力先は下のoutで確認する。
invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
args<-commandArgs(trailingOnly=TRUE);stopifnot(length(args)==1)
base<-normalizePath(args[[1]],winslash="/",mustWork=TRUE)
d<-read.csv(file.path(base,"aggregate_results/matched_key_model_sensitivity.csv"),check.names=FALSE)
# 数据里的 AD 标签在图中显示为 Dementia，比较对象没有改变。
pair_names<-c(CN_vs_AD="CN vs Dementia",CN_vs_MCI="CN vs MCI",MCI_vs_AD="MCI vs Dementia")
biom_names<-c(PTAU_on_ChPICV="P-Tau",sTREM2_on_ChPICV="sTREM2",TAU_on_ChPICV="Tau")
# 各診断ペア内で、同じバイオマーカーのmatchedとunmatchedを並べる。
d<-d[order(match(d$pair_name,names(pair_names)),match(d$model_name,names(biom_names)),match(d$sample_set,c("matched_pair","unmatched_pair"))),]
d$label<-paste(pair_names[d$pair_name],biom_names[d$model_name],ifelse(d$sample_set=="matched_pair","matched","unmatched"),sep=" | ")
out<-file.path(base,"英文投稿包/Supplementary_figures")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
draw<-function(scale=1) {
  n<-nrow(d);xr<-range(c(d$std_low,d$std_high,0));pad<-diff(xr)*.08;xr<-xr+c(-pad,pad)
  par(family="sans",mar=c(13.2,76,1.8,3),xaxs="i",yaxs="i",mgp=c(9.2,4.1,0),col="black",fg="black")
  plot(NA,xlim=xr,ylim=c(.5,n+.5),axes=FALSE,xlab="Standardized beta",ylab="",cex.lab=3.95)
  ticks<-pretty(xr,n=5);ticks<-ticks[ticks>=xr[1]&ticks<=xr[2]]
  axis(1,at=ticks,labels=sprintf("%.2f",ticks),cex.axis=3.85,lwd=2.1*scale,lwd.ticks=2.1*scale,tck=-.018)
  axis(2,at=seq_len(n),labels=rev(d$label),las=1,tick=FALSE,cex.axis=3.85*.92,line=.35)
  box(lwd=1.9*scale)
  segments(0,.5,0,n+.5,lty=2,lwd=1.8*scale,col="gray45")
  yy<-rev(seq_len(n))
# 横线为标准化系数的置信区间，圆点为点估计。
  segments(d$std_low,yy,d$std_high,yy,col="#1F5A8A",lwd=10.5*scale)
  points(d$std_beta,yy,pch=19,col="#1F5A8A",cex=7.2)
}
png(file.path(out,"Supplementary_Figure_S11.png"),width=7800,height=7000,res=260,pointsize=12)
draw();invisible(dev.off())
# True vector export with the same layout at the manuscript's 15.494 cm width.
# PDFは投稿時の幅に合わせて文字と線を縮小する。画像のDPIだけを変更する処理ではない。
s<-6.1/(7800/260)
cairo_pdf(file.path(out,"Supplementary_Figure_S11.pdf"),width=6.1,height=6.1*7000/7800,pointsize=12*s,family="Arial")
draw(s);invisible(dev.off())
write.csv(d,file.path(base,"aggregate_results/S11_plot_values.csv"),row.names=FALSE,fileEncoding="UTF-8")
cat("S11 PNG and vector PDF exported.\n")
