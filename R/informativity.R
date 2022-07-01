
#' Auxiliary function to determine if a vector is informative.
#'
#' @description Informativity is determined by our ability to set a threshold, a value comprised between the minimum and maximum value of the distribution and that
#' separates two distinct groups of observations: the ones above threshold and the ones under. If the threshold is above 1, cells couldn't be distinguished
#' based on activity, the module is uninformative.
#'
#'
#' @param activity Numeric vector, corresponding to activity scores across cells, containing values between 0 and 1.
#' @return Double, threshold distinguishing cell populations.
#'
#'
.activity_assignmentThreshold <- function(activity){

    #### Check parameters ####
    stopifnot(is.numeric(activity))


    nCells <- length(activity)

    thr<-NULL

    #### Looking for bimodals with density curve adjustment 1 ####
    if(is.null(thr)){
        densCurve <- density(activity, adjust=1, cut=0)
        inflPoints <- diff(sign(diff(densCurve$y)))
        # get maxs and mins
        maximumsDens <- which(inflPoints==-2)
        max_y <- densCurve$y[maximumsDens]
        max_x <- densCurve$x[maximumsDens]
        minimumDens <- which(inflPoints==2)
        min_y <- densCurve$y[minimumDens]
        min_x <- densCurve$x[minimumDens]
        # consider only what is after the highest max
        new_start_point<-which(max_y==max(max_y))
        max_x<-max_x[new_start_point:length(max_x)]
        max_y<-max_y[new_start_point:length(max_y)]
        min_x<-min_x[new_start_point:length(min_x)]
        min_y<-min_y[new_start_point:length(min_y)]
        # different cases
        if(length(max_x)>=2){
            temp_thr<-c()
            for(i in 1:(length(max_y)-1)){
                m1<-max_x[i]
                m2<-max_x[i+1]
                min<-min(min_y[which(min_x<m2&min_x>m1)],na.rm=T)
                if(abs(m1-m2)>0 & abs(max(c(max_y[i],max_y[i+1]),na.rm=T)-min)>0.1*max(max_y,na.rm=T)){
                    temp_thr<-c(temp_thr,min_x[which(min_y==min)])
                }
            }
            if(!is.null(temp_thr)){
                thr<-max(temp_thr,na.rm = T)
            }
        }
    }

    #### Looking for bimodals with density curve adjustment 0.5 ####
    if(is.null(thr)){
        densCurve <- density(activity, adjust=0.5, cut=0)
        inflPoints <- diff(sign(diff(densCurve$y)))
        # get maxs and mins
        maximumsDens <- which(inflPoints==-2)
        max_y <- densCurve$y[maximumsDens]
        max_x <- densCurve$x[maximumsDens]
        minimumDens <- which(inflPoints==2)
        min_y <- densCurve$y[minimumDens]
        min_x <- densCurve$x[minimumDens]
        # consider only what is after the highest max
        new_start_point<-which(max_y==max(max_y))
        max_x<-max_x[new_start_point:length(max_x)]
        max_y<-max_y[new_start_point:length(max_y)]
        min_x<-min_x[new_start_point:length(min_x)]
        min_y<-min_y[new_start_point:length(min_y)]
        # different cases
        if(length(max_x)>=2){
            temp_thr<-c()
            for(i in 1:(length(max_y)-1)){
                m1<-max_x[i]
                m2<-max_x[i+1]
                min<-min(min_y[which(min_x<m2&min_x>m1)],na.rm=T)
                if(abs(m1-m2)>0 & abs(max(c(max_y[i],max_y[i+1]),na.rm=T)-min)>0.1*max(max_y,na.rm=T)){
                    temp_thr<-c(temp_thr,min_x[which(min_y==min)])
                }
            }
            if(!is.null(temp_thr)){
                thr<-max(temp_thr,na.rm = T)
            }
        }
    }



    ### Default threshold if none could be attributed ###
    if(is.null(thr)){thr<-1.4}

    return(thr)
}

