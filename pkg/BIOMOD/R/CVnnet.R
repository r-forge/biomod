.CV.nnet = function(Input, Target, size=c(2,4,6, 8), decay=c(0.001, 0.01, 0.05, 0.1), nbCV=5, W=NULL){

    Eval = data.frame(matrix(0, ncol=3, nrow=16, dimnames=list(NULL, c("Size", "Decay", "AUC"))))
    Eval[,1] = rep(size,4)
    Eval[,2] = rep(decay, each=4)
    for(i in 1:nbCV){
        set.seed(555)
        Samp = .SampleMat2(Target, 0.5)
        
        if(is.null(W)){
            Eval[,3] = Eval[,3] + apply(Eval[,1:2], 1, Samp, Target, Input, FUN=function(x, Samp, Target, Input){
              nn = nnet(eval(parse(text = paste("Target[Samp$calibration]",
                    paste(.scopeExpSyst(Biomod.material$VarNames, Biomod.material$VarTypes, "GBM"), collapse = "")))),data=Input[Samp$calibration, ],
                    size = x[1], decay = x[2], maxit = 200, trace = FALSE)
              AUC = .somers2(predict(nn, Input[Samp$evaluation,]), Target[Samp$evaluation])["C"]
              return(AUC)
            })
        } else{
            Eval[,3] = Eval[,3] + apply(Eval[,1:2], 1, Samp, Target, Input, W, FUN=function(x, Samp, Target, Input, W){
              nn = nnet(eval(parse(text = paste("Target[Samp$calibration]",
                    paste(.scopeExpSyst(Biomod.material$VarNames, Biomod.material$VarTypes, "GBM"), collapse = "")))),data=Input[Samp$calibration, ],
                    weights=W[Samp$calibration], size = x[1], decay = x[2], maxit = 200, trace = FALSE)
              AUC = .somers2(predict(nn, Input[Samp$evaluation,]), Target[Samp$evaluation])["C"]
              return(AUC)
            })
        }
    }
    Eval[,3] = Eval[,3]/nbCV
    z =which.max(Eval[,3])
    return(Eval[z, 1:2])
}
