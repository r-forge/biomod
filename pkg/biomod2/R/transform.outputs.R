setGeneric( ".transform.outputs", 
            def = function(modOut, out = 'evaluation',...){
                    standardGeneric( ".transform.outputs" )
                    } )

setMethod('.transform.outputs', signature(modOut='array'), 
  function(modOut, out = 'evaluation'){
    # check out attr
    if(!(out %in% c('evaluation', 'prediction', 'var.import', 'calib.failure', 'models.run', 'prediction.eval' ) )){
      stop(paste("out argument must be one of ", toString(c('evaluation', 'prediction', 'var.import',
                                                            'calib.failure', 'models.run', 'prediction.eval' ))))
    }
    
    # check dim of input list
    if(length(dim(modOut)) != 4 ){
      cat('\n',dim(modOut),'\n')
      print(dimnames(modOut))
      warning("Not computed .transform.outputs because of an imcompatible input list dimention", immediate=T)
      return(NULL)
    }

    if(dim(modOut)[4] == 1 & length(unlist(strsplit(unlist(dimnames(modOut)[4]),'_'))) == 1 ){
      dataset.names <- 'AllData'
    } else{
      if(length(dimnames(modOut)[[4]]) > 0){
        dataset.names <- unlist(sapply(unlist(dimnames(modOut)[4]), function(name){return(tail(unlist(strsplit(name,'_')),1))}))
      } else {
        dataset.names <- paste('PA', 1:dim(modOut)[4])
      }
    }

    run.eval.names <- sub('_','',unlist(dimnames(modOut)[3]))
    mod.names <- unlist(dimnames(modOut)[2])
    
    if (out=='evaluation'){
      if( is.null(modOut['evaluation',1,1,1])){ return(NULL) }
      eval.meth.names <- rownames(as.data.frame(modOut['evaluation',1,1,1]))
      eval.col.names <- colnames(as.data.frame(modOut['evaluation',1,1,1]))
  
      eval.out <- array(data = unlist(modOut['evaluation',,,]),
                        dim = c(length(eval.meth.names),
                                length(eval.col.names),
                                length(mod.names),
                                length(run.eval.names),
                                length(dataset.names)),
                        dimnames = list(eval.meth.names,
                                         eval.col.names,
                                         mod.names,
                                         run.eval.names,
                                         dataset.names))
      
     return(eval.out) 
    }
    
    if (out=='prediction'){
      if( is.null(modOut['pred',1,1,1])){ return(NULL) }
      nb.pts.pred <- length(as.numeric(unlist(modOut['pred',1,1,1])))
      pred.out <- array(data = unlist(modOut['pred',,,]),
                        dim = c(nb.pts.pred,
                                length(mod.names),
                                length(run.eval.names),
                                length(dataset.names)),
                        dimnames = list(NULL,
                                         mod.names,
                                         run.eval.names,
                                         dataset.names))
      
     return(pred.out) 
    }
    
    if (out=='prediction.eval'){
      if( is.null(modOut['pred.eval',1,1,1])){ return(NULL) }
      nb.pts.pred.eval <- length(as.numeric(unlist(modOut['pred.eval',1,1,1])))
      pred.eval.out <- array(data = unlist(modOut['pred.eval',,,]),
                        dim = c(nb.pts.pred.eval,
                                length(mod.names),
                                length(run.eval.names),
                                length(dataset.names)),
                        dimnames = list(NULL,
                                         mod.names,
                                         run.eval.names,
                                         dataset.names))
      
     return(pred.eval.out) 
    }
    
    if (out=='var.import'){
      if( is.null(unlist(modOut['var.import',1,1,1]))){ return(NULL) }
      nb.var <- length(as.numeric(unlist(modOut['var.import',1,1,1])))
  
      vi.out <- array(data = unlist(modOut['var.import',,,]),
                        dim = c(nb.var,
                                length(mod.names),
                                length(run.eval.names),
                                length(dataset.names)),
                        dimnames = list(paste('Var',1:nb.var,sep=''), # to change
                                         mod.names,
                                         run.eval.names,
                                         dataset.names))
      
     return(vi.out) 
    }
    
    if (out == 'calib.failure'){
      cf.out <- unlist(modOut['calib.failure',,,])
      return(cf.out[!is.null(cf.out)])
    }
    
    if (out == 'models.run'){
      mod.run.out <- unlist(modOut['ModelName',,,])
      return(mod.run.out[!is.null(mod.run.out)])
    }
    
  })

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
setMethod('.transform.outputs', signature(modOut='list'), 
  function(modOut, out = 'evaluation', dim.names = NULL){
    # check out attr
    if(!(out %in% c('evaluation', 'prediction', 'prediction.eval', 'var.import', 'calib.failure',
                    'models.run', 'EF.prediction', 'EF.PCA.median', 'EF.evaluation' ) )){
      stop(paste("out argument must be one of ", toString(c('evaluation', 'prediction', 'prediction.eval', 'var.import',
                                                            'calib.failure', 'models.run', 'EF.prediction',
                                                            'EF.PCA.median', 'EF.evaluation',  ))))
    }
    
    if(length(modOut) == 1 & length(unlist(strsplit(unlist(names(modOut)),'_'))) == 1 ){
      dataset.names <- 'AllData'
    } else{
      if(is.null(dim.names)){
        dataset.names <- unlist(sapply(unlist(names(modOut)), function(name){return(tail(unlist(strsplit(name,'_')),1))}))
      } else{
        dataset.names <- unlist(dim.names[1])
      }
    }

    if(is.null(dim.names)){
      run.eval.names <- sub('_','',unlist(names(modOut[[1]]))) # may be good here to test that all names are identics
    
      mod.names <- unlist(names(modOut[[1]][[1]]))
    } else{
      run.eval.names <- unlist(dim.names[2])
      mod.names <- unlist(dim.names[3])
    }
    
    if (out=='evaluation'){
      if( is.null(modOut[[1]][[1]][[1]][['evaluation']])){ return(NULL) }
      eval.meth.names <- rownames(as.data.frame(modOut[[1]][[1]][[1]][['evaluation']]))
      eval.col.names <- colnames(as.data.frame(modOut[[1]][[1]][[1]][['evaluation']]))
      
      eval.out <- lapply(names(modOut),function(d1){ # data set
                    lapply(names(modOut[[d1]]), function(d2){ # run eval
                      lapply(names(modOut[[d1]][[d2]]), function(d3){ # models
                        return(data.frame(modOut[[d1]][[d2]][[d3]][['evaluation']]))
                      })
                    })
                  })
  
      eval.out <- array(data = unlist(eval.out),
                        dim = c(length(eval.meth.names),
                                length(eval.col.names),
                                length(mod.names),
                                length(run.eval.names),
                                length(dataset.names)),
                        dimnames = list(eval.meth.names,
                                         eval.col.names,
                                         mod.names,
                                         run.eval.names,
                                         dataset.names))
      
     return(eval.out) 
    }
    
    if (out=='prediction'){
      if( is.null(modOut[[1]][[1]][[1]][['pred']])){ return(NULL) }

      nb.pts.pred <- length(as.numeric(unlist(modOut[[1]][[1]][[1]][['pred']])))
      
      pred.out <- lapply(names(modOut),function(d1){ # data set
                    lapply(names(modOut[[d1]]), function(d2){ # run eval
                      lapply(names(modOut[[d1]][[d2]]), function(d3){ # models
                        return(as.numeric(modOut[[d1]][[d2]][[d3]][['pred']]))
                      })
                    })
                  })
      
      pred.out <- array(data = unlist(pred.out),
                        dim = c(nb.pts.pred,
                                length(mod.names),
                                length(run.eval.names),
                                length(dataset.names)),
                        dimnames = list(NULL,
                                         mod.names,
                                         run.eval.names,
                                         dataset.names))
      
     return(pred.out) 
    }
    
    if (out=='prediction.eval'){
      if( is.null(modOut[[1]][[1]][[1]][['pred.eval']])){ return(NULL) }

      nb.pts.pred.eval <- length(as.numeric(unlist(modOut[[1]][[1]][[1]][['pred.eval']])))
      
      pred.eval.out <- lapply(names(modOut),function(d1){ # data set
                    lapply(names(modOut[[d1]]), function(d2){ # run eval
                      lapply(names(modOut[[d1]][[d2]]), function(d3){ # models
                        return(as.numeric(modOut[[d1]][[d2]][[d3]][['pred.eval']]))
                      })
                    })
                  })
      
      pred.eval.out <- array(data = unlist(pred.eval.out),
                        dim = c(nb.pts.pred.eval,
                                length(mod.names),
                                length(run.eval.names),
                                length(dataset.names)),
                        dimnames = list(NULL,
                                         mod.names,
                                         run.eval.names,
                                         dataset.names))
      
     return(pred.eval.out) 
    }
    
    if (out=='var.import'){
      if( is.null(modOut[[1]][[1]][[1]][['var.import']])){ return(NULL) }
      nb.var <- length(as.numeric(unlist(modOut[[1]][[1]][[1]][['var.import']])))
      
      ef.mod <- grep(pattern="EF.",mod.names) # EF models
      if(length(ef.mod)>0){
        kept.mod <- mod.names[-ef.mod]
      } else{
        kept.mod <- mod.names
      }
                                        
      vi.out <- lapply(names(modOut),function(d1){ # data set
                  lapply(names(modOut[[d1]]), function(d2){ # run eval
                    lapply(kept.mod, function(d3){ # models without EF ones
                      return(as.matrix(modOut[[d1]][[d2]][[d3]][['var.import']]))
                    })
                  })
                })                                 
  
      vi.out <- array(data = unlist(vi.out),
                        dim = c(nb.var,
                                length(kept.mod),
                                length(run.eval.names),
                                length(dataset.names)),
                        dimnames = list(paste('Var',1:nb.var,sep=''), # to change
                                         kept.mod,
                                         run.eval.names,
                                         dataset.names))
      
     return(vi.out) 
    }
    
    if (out == 'calib.failure'){
      cf.out <- lapply(names(modOut),function(d1){ # data set
                  lapply(names(modOut[[d1]]), function(d2){ # run eval
                    lapply(names(modOut[[d1]][[d2]]), function(d3){ # models 
                      return(as.numeric(modOut[[d1]][[d2]][[d3]][['calib.failure']]))
                    })
                  })
                })
      cf.out <- na.omit(unlist(cf.out))
      cf.out <- cf.out[!is.null(cf.out)]
      if(is.na(cf.out)) cf.out <- 'none'
      return(cf.out)
    }
    
    if (out == 'models.run'){
      mod.run.out <- lapply(names(modOut),function(d1){ # data set
                  lapply(names(modOut[[d1]]), function(d2){ # run eval
                    lapply(names(modOut[[d1]][[d2]]), function(d3){ # models 
                      return(as.character(modOut[[d1]][[d2]][[d3]][['ModelName']]))
                    })
                  })
                })  
      mod.run.out <- na.omit(unlist(mod.run.out))
      return(mod.run.out[!is.null(mod.run.out)])
    }
    

    if (out == 'EF.prediction'){
      if( is.null(modOut[[1]][[1]][[1]][['EM']])){ return(NULL) }

      nb.pts.ef.pred <- length(as.numeric(unlist(modOut[[1]][[1]][[1]][['EM']])))
      
      ef.pred.out <- lapply(1:length(modOut),function(d1){ # data set
                        lapply(1:length(modOut[[d1]]), function(d2){ # run eval
                          lapply(1:length(modOut[[d1]][[d2]]), function(d3){ # models
                            return(as.numeric(modOut[[d1]][[d2]][[d3]][['EM']]))
                          })
                        })
                      })
      
      ef.pred.out <- array( data = unlist(ef.pred.out),
                            dim = c(nb.pts.ef.pred,
                                    length(modOut[[1]][[1]]),
                                    length(modOut[[1]]),
                                    length(modOut)),
                            dimnames = list(NULL,
                                             mod.names,
                                             run.eval.names,
                                             dataset.names))
      
     return(ef.pred.out) 
    }

    if (out == 'EF.PCA.median'){
      if( is.null(modOut[[1]][[1]][[1]][['PCA.median']])){ return(NULL) }

      ef.pca.out <- lapply(1:length(modOut),function(d1){ # data set
                        lapply(1:length(modOut[[d1]]), function(d2){ # run eval
                          lapply(1:length(modOut[[d1]][[d2]]), function(d3){ # models
                            return(as.character(modOut[[d1]][[d2]][[d3]][['PCA.median']]))
                          })
                        })
                      })
      
      ef.pca.out <- array( data = unlist(ef.pca.out),
                            dim = c(1,
                                    length(modOut[[1]][[1]]),
                                    length(modOut[[1]]),
                                    length(modOut)),
                            dimnames = list(NULL,
                                             mod.names,
                                             run.eval.names,
                                             dataset.names))
      
     return(ef.pca.out) 
    }

    if (out == 'EF.evaluation'){
      if( is.null(modOut[[1]][[1]][[1]][['EM.eval']])){ return(NULL) }
      eval.meth.names <- rownames(as.data.frame(modOut[[1]][[1]][[1]][['EM.eval']]))
      eval.col.names <- colnames(as.data.frame(modOut[[1]][[1]][[1]][['EM.eval']]))
      
      ef.eval.out <- lapply(1:length(modOut),function(d1){ # data set
                    lapply(1:length(modOut[[d1]]), function(d2){ # run eval
                      lapply(1:length(modOut[[d1]][[d2]]), function(d3){ # models
                        return(data.frame(modOut[[d1]][[d2]][[d3]][['EM.eval']]))
                      })
                    })
                  })
  
      ef.eval.out <- array(data = unlist(ef.eval.out),
                        dim = c(length(eval.meth.names),
                                length(eval.col.names),
                                length(modOut[[1]][[1]]),
                                length(modOut[[1]]),
                                length(modOut)),
                        dimnames = list(eval.meth.names,
                                         eval.col.names,
                                         mod.names,
                                         run.eval.names,
                                         dataset.names))
      
     return(ef.eval.out)
    }
    
  })

DF_to_ARRAY <- function(df){
  if(!is.data.frame(df) & !is.matrix(df)){
    if(is.list(df)){
      df.names <- names(df)
      df <- as.data.frame(df)
      names(df) <- df.names
    } else{
      stop("You have to give a data.frame")
    }
  }
  
  array.dim.names <- c(list(c(NULL)),rev(apply(sapply(strsplit(colnames(df), '_'), tail, n=3),1,unique)))
  array.dim <- c(nrow(df),sapply(array.dim.names[-1],length))
  array.out <- array(data=NA, dim=array.dim, dimnames=array.dim.names)

  for(x in colnames(df)){
    dimTmp <- rev(tail(unlist(strsplit(x, '_')), n=3))
    array.out[,dimTmp[1],dimTmp[2],dimTmp[3]] <- df[,x]
  }
  return(array.out)
}