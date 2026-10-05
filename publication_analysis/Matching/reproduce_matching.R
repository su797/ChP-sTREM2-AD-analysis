# Usage: Rscript reproduce_matching.R DATA_CLEAN_DIR TRANSFORMATION_PLAN_CSV OUTPUT_DIR
# Inputs remain local; output contains aggregate estimates only.
# 此脚本读取已经生成的匹配样本，不在这里重新匹配。
if (.Platform$OS.type == "windows") {
  invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
}
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=3L)stop("Supply DATA_CLEAN_DIR TRANSFORMATION_PLAN_CSV OUTPUT_DIR")
script_arg<-grep("^--file=",commandArgs(FALSE),value=TRUE)[1]
script_dir<-dirname(normalizePath(sub("^--file=","",script_arg),winslash="/"))
source(file.path(script_dir,"s11_model_functions.R"),encoding="UTF-8")
project_config<-dget(file.path(script_dir,"s11_model_config.R"))
result_data_clean_dir<-normalizePath(args[1],winslash="/",mustWork=TRUE)
analysis_data<-read_project_data(file.path(result_data_clean_dir,"ChpSTREM2AD_analysis_dataset.csv"))
# 変換表は元の解析と同じものを指定する。列名だけでなく、各変数の解析尺度もこれで決まる。
transformation_plan<-read_project_data(args[2])
rows<-run_matched_key_models()
# 置信区间和系数使用同一个标准化比例，端点再用 pmin/pmax 排序。
ratio<-rows$std_beta/rows$estimate
rows$std_low<-pmin(rows$conf.low*ratio,rows$conf.high*ratio)
rows$std_high<-pmax(rows$conf.low*ratio,rows$conf.high*ratio)
dir.create(args[3],recursive=TRUE,showWarnings=FALSE)
write.csv(rows,file.path(args[3],"matched_key_model_sensitivity.csv"),row.names=FALSE,fileEncoding="UTF-8")
writeLines(capture.output(sessionInfo()),file.path(args[3],"R_session_info.txt"))
cat(nrow(rows),"aggregate model estimates written.\n")
