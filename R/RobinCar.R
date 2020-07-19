#' Main function for RobinCar
#'
#' @param Fun.data A dataframe formulated as data.frame(y, I, strata_cross, x1, x2,..., xp). Specifically, y is the observed outcome;
#' I is the treatment indicator equals 1,...,k; strata_cross is a factor corresponding to the joint levels of Z; x1,...,xp are additional covariates
#' to be adjusted. Any categorical variable in x1,...,xp needs to be dummied as 0 or 1 (numeric).
#' @param Fun.trt_label.pair A vector of length two, e.g. Fun.trt_label.pair=c(s,t), indicates which two treatments are being compared.
#' The estimand is therefore the average treatment effect comparing the latter to the former.
#' @param Fun.trt_alc.pair A vector of length two, e.g., Fun.trt_alc.pair=c(pi_s, pi_t), that corresponds to the target probability of taking treatments s and t, respectively.
#' If Fun.trt_alc.pair is NA, pi_s and pi_t will be estimated from Fun.data$I.
#' @param Fun.covariates A subset of c('x1', ..., 'xp'), which indicates the covariates used for adjustment and should match the column names in Fun.data.
#' The default is NA, that is no covariate adjustment is made and theta_hat will be applied.
#' @param Fun.adjustment Either "Hetero" or "Homo". "Hetero" corresponds to theta_A_hat; "Homo" corresponds to theta_B_hat.
#' If Fun.covariates is NA, this argument will be omitted. For recommendations of using either "Hetero" or "Homo", please refer to the last section of Ye et al., (2020).
#'
#' @return A list
#' \describe{
#' \item{Z_statistic}{Equals \code{estimate}/\code{se}}
#' \item{estimate}{Estimated average treatment effect}
#' \item{se}{Standard error of \code{estimate}}
#' \item{Pval}{Two-sided p-value calculated using the \code{Z_statistic}}
#' }
#' @references Ting Ye, Yanyao Yi, Jun Shao (2020). Inference on Average Treatment Effect under Minimization and Other Covariate-Adaptive Randomization Methods.
#'
#' @import stats tidyverse
#' @export
#'
RobinCar<-function(Fun.data, Fun.trt_label.pair, Fun.trt_alc.pair = NA,
                    Fun.covariates = NA, Fun.adjmethod="Hetero"){
  # print("For now only two sided pval is provided, one sided pval can be calculated from results!")
  if(! is.data.frame(Fun.data)){
    stop(print("Fun.data should be a data frame!"))
  }
  if(! "strata_cross" %in% names(Fun.data)){
    stop(print("Please name the strata column of Fun.data as 'strata_cross'!"))
  }else if(! is.factor(Fun.data$strata_cross)){
    stop(print("strata_cross should be a factor, even if there is only one level, i.e., no stratification factor!"))
  }
  if(! "y" %in% names(Fun.data)){
    stop(print("Please name the outcome column of Fun.data as 'y'!"))
  }
  if(! "I" %in% names(Fun.data)){
    stop(print("Please name the treatment assignment column of Fun.data as 'I'!"))
  }
  if(is.factor(Fun.data$I)){
    if(! all(Fun.trt_label.pair %in% levels(Fun.data$I))){
      stop(print("Please make sure Fun.data$I contains Fun.trt_label.pair!"))
    }
  }else{
    if(! all(Fun.trt_label.pair %in% Fun.data$I)){
      stop(print("Please make sure Fun.data$I contains Fun.trt_label.pair!"))
    }
  }
  if(length(Fun.trt_label.pair)!=2){
    stop(print("Sorry, we can only do pairwise comparison once a time for now! Please make sure Fun.trt_label.pair is of length two."))
  }
  if(is.na(Fun.covariates)){
    print("No covariate adjustment is applied!")
  }else if(Fun.covariates %in% colnames(Fun.data)){
    print(paste0(Fun.adjmethod, "geneous treatment effect for ", Fun.covariates, " is considered!"))
  }else{
    stop(print("Fun.covariates does match colnames of Fun.data! Please check and match the name of covariates!"))
  }
  if(min(table(Fun.data[, c("I", "strata_cross")])) < 10){
    warning("The smallest sample size within a strata and treatment is fewer than 10!")
  }


  Fun.n <- nrow(Fun.data)

  if(all(is.na(Fun.trt_alc.pair))){
    Fun.trt_pie <- Fun.data %>% group_by(I) %>% summarise(pieI=n()/Fun.n) %>%
      filter(I %in% Fun.trt_label.pair) %>% arrange(match(Fun.trt_label.pair, I))
  }else{
    Fun.trt_pie <- Fun.data %>% distinct(I) %>%
      filter(I %in% Fun.trt_label.pair) %>% arrange(match(Fun.trt_label.pair, I)) %>%
      bind_cols(pieI = Fun.trt_alc.pair)
  }

  Fun.estimate <- 0
  Fun.sigmaV_sq <- 0
  Fun.sigmaU_sq <- 0
  Fun.sigmaA_sq <- 0
  Fun.sigmaB_sq <- 0


  if(all(is.na(Fun.covariates))){
    # no adjustment
    for(Fun.j in levels(Fun.data$strata_cross)){
      Fun.data.z <- Fun.data[Fun.data$strata_cross==Fun.j, ]
      Fun.nz <- nrow(Fun.data.z)
      Fun.sumstat.z <- Fun.data.z %>%
        group_by(I) %>%
        summarise(ybar = mean(y), yvar = var(y)) %>%
        filter(I %in% Fun.trt_label.pair) %>%
        arrange(match(Fun.trt_label.pair, I)) %>%
        left_join(y=Fun.trt_pie, by = "I")
      Fun.estimate <- Fun.estimate + (Fun.nz/Fun.n)*(Fun.sumstat.z$ybar[2]-Fun.sumstat.z$ybar[1])
      Fun.sigmaU_sq <- Fun.sigmaU_sq + (Fun.nz/Fun.n)*sum(Fun.sumstat.z$yvar/Fun.sumstat.z$pieI)
      Fun.sigmaV_sq <- Fun.sigmaV_sq + (Fun.nz/Fun.n)*(Fun.sumstat.z$ybar[2]-Fun.sumstat.z$ybar[1])^2
    }
    Fun.sd <- sqrt(Fun.sigmaU_sq + Fun.sigmaV_sq - Fun.estimate^2)

  }else if(Fun.adjmethod=="Hetero"){
    # method A
    for(Fun.j in levels(Fun.data$strata_cross)){
      Fun.data.z <- Fun.data[Fun.data$strata_cross==Fun.j, ]
      Fun.nz <- nrow(Fun.data.z)
      Fun.data.z.model <- Fun.data.z %>%
        mutate_at(Fun.covariates, list(centered = ~ scale(., scale = F))) %>%
        select(y, I, contains("centered")) %>%
        filter(I %in% Fun.trt_label.pair) %>%
        mutate(I_st = as.numeric(factor(I, levels = Fun.trt_label.pair))-1) %>%
        mutate_at(vars(contains("centered")), list(beta1 = ~. * I_st, beta0 = ~. * (1-I_st))) %>%
        select(y, I_st, contains("beta"))
      Fun.model.z.lm <- lm(y~., data = Fun.data.z.model)
      Fun.estimate <- Fun.estimate + (Fun.nz/Fun.n)*Fun.model.z.lm$coefficients[2]
      Fun.beta1 <- Fun.model.z.lm$coefficients[str_detect(names(Fun.model.z.lm$coefficients), "beta1")]
      Fun.beta0 <- Fun.model.z.lm$coefficients[str_detect(names(Fun.model.z.lm$coefficients), "beta0")]
      Fun.sigmaA_sq_2 <- (Fun.beta1-Fun.beta0)%*%var(Fun.data.z[, Fun.covariates])%*%(Fun.beta1-Fun.beta0)
      Fun.sigmaA_sq_1 <- sum(c(var(Fun.data.z$y[Fun.data.z$I==Fun.trt_label.pair[1]]
                                     - rowSums(sweep(as.matrix(Fun.data.z[Fun.data.z$I==Fun.trt_label.pair[1], Fun.covariates]), MARGIN = 2, Fun.beta0, "*"))),
                                 var(Fun.data.z$y[Fun.data.z$I==Fun.trt_label.pair[2]]
                                     - rowSums(sweep(as.matrix(Fun.data.z[Fun.data.z$I==Fun.trt_label.pair[2], Fun.covariates]), MARGIN = 2, Fun.beta1, "*")))
      )/Fun.trt_pie$pieI)
      Fun.sigmaA_sq <- Fun.sigmaA_sq + (Fun.nz/Fun.n)*as.vector(Fun.sigmaA_sq_1+Fun.sigmaA_sq_2)
      Fun.sigmaV_sq <- Fun.sigmaV_sq + (Fun.nz/Fun.n)*as.vector(Fun.model.z.lm$coefficients[2])^2
    }
    Fun.sd <- sqrt(Fun.sigmaA_sq + Fun.sigmaV_sq - Fun.estimate^2)

  }else if(Fun.adjmethod=="Homo"){
    # method B
    for(Fun.j in levels(Fun.data$strata_cross)){
      Fun.data.z <- Fun.data[Fun.data$strata_cross==Fun.j, ]
      Fun.nz <- nrow(Fun.data.z)
      Fun.data.z.model <- Fun.data.z %>%
        mutate(I_st = factor(I)) %>%
        select(y, I_st, Fun.covariates)
      Fun.model.z.lm <- lm(y~0+., data = Fun.data.z.model)
      Fun.thetahat.Bz <- diff(Fun.model.z.lm$coefficients[match(Fun.trt_label.pair, levels(Fun.data.z.model$I_st))])
      Fun.estimate <- Fun.estimate + (Fun.nz/Fun.n)*Fun.thetahat.Bz
      Fun.betahat <- Fun.model.z.lm$coefficients[-(1:length(levels(Fun.data.z.model$I_st)))]
      Fun.sigmaB_sq <- Fun.sigmaB_sq +
        sum(c(var(Fun.data.z$y[Fun.data.z$I==Fun.trt_label.pair[1]]
                  - rowSums(sweep(as.matrix(Fun.data.z[Fun.data.z$I==Fun.trt_label.pair[1], Fun.covariates]), MARGIN = 2, Fun.betahat, "*"))),
              var(Fun.data.z$y[Fun.data.z$I==Fun.trt_label.pair[2]]
                  - rowSums(sweep(as.matrix(Fun.data.z[Fun.data.z$I==Fun.trt_label.pair[2], Fun.covariates]), MARGIN = 2, Fun.betahat, "*")))
        )/Fun.trt_pie$pieI)*(Fun.nz/Fun.n)
      Fun.sigmaV_sq <- Fun.sigmaV_sq + (Fun.nz/Fun.n)*as.vector(Fun.thetahat.Bz)^2
    }
    Fun.sd <- sqrt(Fun.sigmaB_sq + Fun.sigmaV_sq - Fun.estimate^2)

  }else{
    stop('There is no such covariate adjustment method, please choose from "Hetero" or "Homo"!')
  }


  return(list(Z_statistic = as.vector(sqrt(Fun.n)*Fun.estimate/Fun.sd),
              estimate = as.vector(Fun.estimate),
              se = as.vector(Fun.sd/sqrt(Fun.n)),
              Pval = as.vector(2*pnorm(q=abs(sqrt(Fun.n)*Fun.estimate/Fun.sd), lower.tail = F))))
}










