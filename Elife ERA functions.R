#----------------------------------------------------------------------------------------------------
# title: "Combination of inflammatory and vascular markers in the febrile phase of dengue is associated with more severe outcomes"
# author: "Nguyen Lam Vuong"
# date: "29-Jul-2021"
# SOME FUNCTIONS FOR USE IN THE ANALYSIS

#----------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------
# Scatter matrix plot
## This is modified from Pascal GP Martin's codes. Thanks him a lot!
## https://pascal-martin.netlify.app/post/nicer-scatterplot-in-gggally/

GGscatterPlot <- function(data, mapping, ..., 
                          method = "spearman") {
  
  #Assemble data frame
  x <- GGally::eval_data_col(data, mapping$x)
  y <- GGally::eval_data_col(data, mapping$y)
  df <- data.frame(x = x, y = y)
  df <- df[!is.na(x) & !is.na(y),]
  #Get correlation coefficient
  cor <- cor(df$x, df$y, method = method)
  # PCA
  nonNull <- x!=0 & y!=0
  dfpc <- prcomp(~x+y, df[nonNull,])
  df$cols <- predict(dfpc, df)[,1]
  # Define the direction of color range based on PC1 orientation:
  dfsum <- x+y
  colDirection <- ifelse(dfsum[which.max(df$cols)] < 
                           dfsum[which.min(df$cols)],
                         1,
                         -1)
  #Get 2D density for alpha
  dens2D <- MASS::kde2d(df$x, df$y)
  df$density <- fields::interp.surface(dens2D , 
                                       df[,c("x", "y")])
  
  if (any(df$density==0)) {
    mini2D = min(df$density[df$density!=0]) #smallest non zero value
    df$density[df$density==0] <- mini2D
  }
  #Prepare plot
  pp <- ggplot(df, aes(x=x, y=y, color = cols, alpha = 1/density)) +
    ggplot2::geom_point(shape=16, size=1, show.legend = FALSE) +
    ggplot2::scale_alpha(range = c(.2, .3)) +
    ggplot2::geom_label(
      data = data.frame(
        xlabel = max(x, na.rm = TRUE),
        ylabel = min(y, na.rm = TRUE),
        lab = round(cor, digits = 2)),
      mapping = ggplot2::aes(x = xlabel, 
                             y = ylabel, 
                             label = lab),
      hjust = 1, vjust = 0, alpha = .2,
      size = 3, fontface = "bold", label.size = NA,
      inherit.aes = F # do not inherit anything from the ...
    ) +
    theme_minimal() +
    theme(axis.text.x=element_text(angle=60))
  return(pp)
}


#----------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------
# Function to customize facet_wrap

scale_override <- function(which, scale) {
  if(!is.numeric(which) || (length(which) != 1) || (which %% 1 != 0)) {
    stop("which must be an integer of length 1")
  }
  
  if(is.null(scale$aesthetics) || !any(c("x", "y") %in% scale$aesthetics)) {
    stop("scale must be an x or y position scale")
  }
  
  structure(list(which = which, scale = scale), class = "scale_override")
}

CustomFacetWrap <- ggproto(
  "CustomFacetWrap", FacetWrap,
  init_scales = function(self, layout, x_scale = NULL, y_scale = NULL, params) {
    # make the initial x, y scales list
    scales <- ggproto_parent(FacetWrap, self)$init_scales(layout, x_scale, y_scale, params)
    
    if(is.null(params$scale_overrides)) return(scales)
    
    max_scale_x <- length(scales$x)
    max_scale_y <- length(scales$y)
    
    # ... do some modification of the scales$x and scales$y here based on params$scale_overrides
    for(scale_override in params$scale_overrides) {
      which <- scale_override$which
      scale <- scale_override$scale
      
      if("x" %in% scale$aesthetics) {
        if(!is.null(scales$x)) {
          if(which < 0 || which > max_scale_x) stop("Invalid index of x scale: ", which)
          scales$x[[which]] <- scale$clone()
        }
      } else if("y" %in% scale$aesthetics) {
        if(!is.null(scales$y)) {
          if(which < 0 || which > max_scale_y) stop("Invalid index of y scale: ", which)
          scales$y[[which]] <- scale$clone()
        }
      } else {
        stop("Invalid scale")
      }
    }
    
    # return scales
    scales
  }
)

facet_wrap_custom <- function(..., scale_overrides = NULL) {
  # take advantage of the sanitizing that happens in facet_wrap
  facet_super <- facet_wrap(...)
  
  # sanitize scale overrides
  if(inherits(scale_overrides, "scale_override")) {
    scale_overrides <- list(scale_overrides)
  } else if(!is.list(scale_overrides) || 
            !all(vapply(scale_overrides, inherits, "scale_override", FUN.VALUE = logical(1)))) {
    stop("scale_overrides must be a scale_override object or a list of scale_override objects")
  }
  
  facet_super$params$scale_overrides <- scale_overrides
  
  ggproto(NULL, CustomFacetWrap,
          shrink = facet_super$shrink,
          params = facet_super$params
  )
}


#----------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------
# Function for splines

k <- function(x) {out <- quantile(x, .5)}
B <- function(x) {out <- c(quantile(x, .1), quantile(x, .9))}
kB <- function(x) {out <- c(quantile(x, .1), quantile(x, .5), quantile(x, .9))}

ns1 <- function(x) {
  require(splines)
  out <- ns(x, k=k(x), B=B(x))
}


#----------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------
# Function for screening (Appendix 3)
screen <- function(mod, bio, dat) {
  if (mod==1) {f <- paste("sev.or.inte ~ ", bio, sep="")}
  else {if (mod==2) {f <- paste("sev.or.inte ~ ns1(", bio, ")", sep="")}
    else {if (mod==3) {f <- paste("sev.or.inte ~ ", bio, " + u", bio, sep="")}
      else {if (mod==4) {f <- paste("sev.or.inte ~ ns1(", bio, ") + u", bio, sep="")}
        else {if (mod==5) {f <- paste("sev.or.inte ~ ", bio, " * Age", sep="")}
          else {if (mod==6) {f <- paste("sev.or.inte ~ ns1(", bio, ") * Age", sep="")}
            else {if (mod==7) {f <- paste("sev.or.inte ~ ", bio, " * Age + u", bio, sep="")}
              else {f <- paste("sev.or.inte ~ ns1(", bio, ") * Age + u", bio, sep="")}}}}}}}
  if (bio %in% c("SDC", "IL1RA", "Fer", "CRP") & mod %in% c(3,4,7,8)) {out <- NA} 
  else {m <- glm(formula(f), data=dat, family=binomial); out <- round(AIC(m),1)}
  return(out)
}


#----------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------
# Function to get results from models

#----------------------------------------------------------------------------------------------------
## Single model

get_est1 <- function(out, bio, ref1, ref2, est, age, dat) {

  # Perform model
  if (bio %in% c("VCAM", "Ang", "IL8", "IP10", "TREM", "Vir")) {
    if (out=="sev.only") {
      f <- paste(out, " ~ ", bio, " * Age + u", bio, sep="")
    } else {
      f <- paste(out, " ~ rcs(", bio, ", kB(", bio, ")) * rcs(Age, kB(Age)) + u", bio, sep="")
    }
  } else {
    if (out=="sev.only") {
      f <- paste(out, " ~ ", bio, " * Age", sep="")
    } else {
      f <- paste(out, " ~ rcs(", bio, ", kB(", bio, ")) * rcs(Age, kB(Age))", sep="")
    }
  }
  
  if (out=="sev.or.inte") {
    m <- lrm(formula(f), data=dat)
  } else {
    m <- robcov(lrm(formula(f), data=dat, weights=ipwsd, x=T, y=T))
  }
  
  # Get p-values for (1) main effect, (2) interaction, and (3) non-linear effect
  if (est=="p") {output <- anova(m)[1,3]}
  else {
    if (est=="p int") {output <- anova(m)[2,3]}
    else {
      if (est=="p ns") {output <- anova(m)[3,3]}
      
      # Get OR and 95% CI
      else {
        text <- paste("summary(m, ", bio, "=c(", ref1, ",", ref2, "), ", "Age=", age, ")", sep="")
        s <- eval(parse(text = text))
        if (est=="OR") {output <- s[2,4]}
        else {
          if (est=="loCI") {output <- s[2,6]}
          else {output <- s[2,7]}
        }
      }
    }
  }
  return(output)
}

#----------------------------------------------------------------------------------------------------
## Global model without viremia

get_est2 <- function(out, bio, ref1, ref2, est, age, dat) {

  # Perform model
  if (out=="sev.only") {
    f <- paste(out, "~ (VCAM + SDC + Ang + IL8 + IP10 + IL1RA + CD163 + TREM + Fer + CRP) * Age + uVCAM + uAng + uIL8 + uIP10 + uTREM")
  } else {
    f <- paste(out, "~ (rcs(VCAM, parms=kB(VCAM)) + rcs(SDC, parms=kB(SDC)) + rcs(Ang, parms=kB(Ang)) + rcs(IL8, parms=kB(IL8)) + rcs(IP10, parms=kB(IP10)) + rcs(IL1RA, parms=kB(IL1RA)) + rcs(CD163, parms=kB(CD163)) + rcs(TREM, parms=kB(TREM)) + rcs(Fer, parms=kB(Fer)) + rcs(CRP, parms=kB(CRP))) * rcs(Age, kB(Age)) + uVCAM + uAng + uIL8 + uIP10 + uTREM")
  }
  
  if (out=="sev.or.inte") {
    m <- lrm(formula(f), data=dat)
  } else {
    m <- robcov(lrm(formula(f), data=dat, weights=ipwsd, x=T, y=T))
  }
  
  mp <- anova(m)
  
  # Get p-values for (1) main effect, (2) interaction, and (3) non-linear effect
  if (est=="p") {output <- mp[which(rownames(mp)==paste(bio, "(Factor+Higher Order Factors)", sep="  ")), 3]}
  else {
    if (est=="p int") {output <- mp[which(rownames(mp)==paste(bio, "(Factor+Higher Order Factors)", sep="  "))+1, 3]}
    else {
      if (est=="p ns") {output <- mp[which(rownames(mp)==paste(bio, "(Factor+Higher Order Factors)", sep="  "))+2, 3]}
      
      # Get OR and 95% CI
      else {
        text <- paste("summary(m, ", bio, "=c(", ref1, ",", ref2, "), ", "Age=", age, ")", sep="")
        s <- eval(parse(text = text))
        s0 <- s[which(rownames(s)==bio)+1,]
        if (est=="OR") {output <- s0[4]}
        else {
          if (est=="loCI") {output <- s0[6]}
          else {output <- s0[7]}
        }
      }
    }
  }
  return(output)
}

#----------------------------------------------------------------------------------------------------
## Global model with viremia

get_est2v <- function(out, bio, ref1, ref2, est, age, dat) {

  # Perform model
  if (out=="sev.only") {
    f <- paste(out, "~ (VCAM + SDC + Ang + IL8 + IP10 + IL1RA + CD163 + TREM + Fer + CRP + Vir) * Age + uVCAM + uAng + uIL8 + uIP10 + uTREM + uVir")
  } else {
    f <- paste(out, "~ (rcs(VCAM, parms=kB(VCAM)) + rcs(SDC, parms=kB(SDC)) + rcs(Ang, parms=kB(Ang)) + rcs(IL8, parms=kB(IL8)) + rcs(IP10, parms=kB(IP10)) + rcs(IL1RA, parms=kB(IL1RA)) + rcs(CD163, parms=kB(CD163)) + rcs(TREM, parms=kB(TREM)) + rcs(Fer, parms=kB(Fer)) + rcs(CRP, parms=kB(CRP)) + rcs(Vir, parms=kB(Vir))) * rcs(Age, kB(Age)) + uVCAM + uAng + uIL8 + uIP10 + uTREM + uVir")
  }
  
  if (out=="sev.or.inte") {
    m <- lrm(formula(f), data=dat)
  } else {
    m <- robcov(lrm(formula(f), data=dat, weights=ipwsd, x=T, y=T))
  }
  
  mp <- anova(m)
  
  # Get p-values for (1) main effect, (2) interaction, and (3) non-linear effect
  if (est=="p") {output <- mp[which(rownames(mp)==paste(bio, "(Factor+Higher Order Factors)", sep="  ")), 3]}
  else {
    if (est=="p int") {output <- mp[which(rownames(mp)==paste(bio, "(Factor+Higher Order Factors)", sep="  "))+1, 3]}
    else {
      if (est=="p ns") {output <- mp[which(rownames(mp)==paste(bio, "(Factor+Higher Order Factors)", sep="  "))+2, 3]}
      
      # Get OR and 95% CI
      else {
        text <- paste("summary(m, ", bio, "=c(", ref1, ",", ref2, "), ", "Age=", age, ")", sep="")
        s <- eval(parse(text = text))
        s0 <- s[which(rownames(s)==bio)+1,]
        if (est %in% c("OR", "HR")) {output <- s0[4]}
        else {
          if (est=="loCI") {output <- s0[6]}
          else {output <- s0[7]}
        }
      }
    }
  }
  return(output)
}


#----------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------
# Function for getting predicted values from models

#----------------------------------------------------------------------------------------------------
## Single model

get_pred1 <- function(out, bio, age, dat) {

  # Choose variable
  if (bio %in% c("VCAM", "Ang", "IL8", "IP10", "TREM", "Vir")) {
    if (out=="sev.only") {
      f <- paste(out, " ~ ", bio, " * Age + u", bio, sep="")
    } else {
      f <- paste(out, " ~ rcs(", bio, ", kB(", bio, ")) * rcs(Age, kB(Age)) + u", bio, sep="")
    }
  } else {
    if (out=="sev.only") {
      f <- paste(out, " ~ ", bio, " * Age", sep="")
    } else {
      f <- paste(out, " ~ rcs(", bio, ", kB(", bio, ")) * rcs(Age, kB(Age))", sep="")
    }
  }
  
  # Perform model
  if (out=="sev.or.inte") {
    m <- lrm(formula(f), data=dat)
  } else {
    m <- robcov(lrm(formula(f), data=dat, weights=ipwsd, x=T, y=T))
  }
  
  # Get prediction
  text <- paste("Predict(m, ", bio, "=dat$", bio, ", Age=", age, ", ref.zero=T)", sep="")
  p <- eval(parse(text = text))
  output <- as.data.frame(p) %>% mutate(biomarker = bio, model = "Single model")
  return(output)
}

#----------------------------------------------------------------------------------------------------
## Global model without viremia

get_pred2 <- function(out, bio, age, dat) {

  # Perform model
  if (out=="sev.only") {
    f <- paste(out, "~ (VCAM + SDC + Ang + IL8 + IP10 + IL1RA + CD163 + TREM + Fer + CRP) * Age + uVCAM + uAng + uIL8 + uIP10 + uTREM")
  } else {
    f <- paste(out, "~ (rcs(VCAM, parms=kB(VCAM)) + rcs(SDC, parms=kB(SDC)) + rcs(Ang, parms=kB(Ang)) + rcs(IL8, parms=kB(IL8)) + rcs(IP10, parms=kB(IP10)) + rcs(IL1RA, parms=kB(IL1RA)) + rcs(CD163, parms=kB(CD163)) + rcs(TREM, parms=kB(TREM)) + rcs(Fer, parms=kB(Fer)) + rcs(CRP, parms=kB(CRP))) * rcs(Age, kB(Age)) + uVCAM + uAng + uIL8 + uIP10 + uTREM")
  }
  
  if (out=="sev.or.inte") {
    m <- lrm(formula(f), data=dat)
  } else {
    m <- robcov(lrm(formula(f), data=dat, weights=ipwsd, x=T, y=T))
  }
  
  # Get prediction
  text <- paste("Predict(m, ", bio, "=dat$", bio, ", Age=", age, ", ref.zero=T)", sep="")
  p <- eval(parse(text = text))
  output <- as.data.frame(p) %>% mutate(biomarker = bio, model = "Global model")
  return(output)
}

#----------------------------------------------------------------------------------------------------
## Global model with viremia

get_pred2v <- function(out, bio, age, dat) {

  # Perform model
  if (out=="sev.only") {
    f <- paste(out, "~ (VCAM + SDC + Ang + IL8 + IP10 + IL1RA + CD163 + TREM + Fer + CRP + Vir) * Age + uVCAM + uAng + uIL8 + uIP10 + uTREM + uVir")
  } else {
    f <- paste(out, "~ (rcs(VCAM, parms=kB(VCAM)) + rcs(SDC, parms=kB(SDC)) + rcs(Ang, parms=kB(Ang)) + rcs(IL8, parms=kB(IL8)) + rcs(IP10, parms=kB(IP10)) + rcs(IL1RA, parms=kB(IL1RA)) + rcs(CD163, parms=kB(CD163)) + rcs(TREM, parms=kB(TREM)) + rcs(Fer, parms=kB(Fer)) + rcs(CRP, parms=kB(CRP)) + rcs(Vir, parms=kB(Vir))) * rcs(Age, kB(Age)) + uVCAM + uAng + uIL8 + uIP10 + uTREM + uVir")
  }
  
  if (out=="sev.or.inte") {
    m <- lrm(formula(f), data=dat)
  } else {
    m <- robcov(lrm(formula(f), data=dat, weights=ipwsd, x=T, y=T))
  }
  
  # Get prediction
  text <- paste("Predict(m, ", bio, "=dat$", bio, ", Age=", age, ", ref.zero=T)", sep="")
  p <- eval(parse(text = text))
  output <- as.data.frame(p) %>% mutate(biomarker = bio, model = "Global model")
  return(output)
}


#----------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------
# Function for getting predicted values from models with the interaction with serotype

#----------------------------------------------------------------------------------------------------
## Single model

get_pred1s <- function(out, bio, age, serotype, dat) {

  # Choose variable
  if (bio %in% c("VCAM", "Ang", "IL8", "IP10", "TREM", "Vir")) {
    if (out=="sev.only") {
      f <- paste(out, " ~ ", bio, " * Age * Serotype + u", bio, sep="")
    } else {
      f <- paste(out, " ~ rcs(", bio, ", kB(", bio, ")) * rcs(Age, kB(Age)) + rcs(", bio, ", kB(", bio, ")) * Serotype + u", bio, sep="")
    }
  } else {
    if (out=="sev.only") {
      f <- paste(out, " ~ ", bio, " * Age", sep="")
    } else {
      f <- paste(out, " ~ rcs(", bio, ", kB(", bio, ")) * rcs(Age, kB(Age)) + rcs(", bio, ", kB(", bio, ")) * Serotype", sep="")
    }
  }
  
  # Perform model
  if (out=="sev.or.inte") {
    m <- lrm(formula(f), data=dat)
  } else {
    m <- robcov(lrm(formula(f), data=dat, weights=ipwsd, x=T, y=T))
  }
  
  # Get prediction
  text <- paste("Predict(m, ", bio, "=dat$", bio, ", Age=", age, ", Serotype='", serotype, "', ref.zero=T)", sep="")
  p <- eval(parse(text = text))
  output <- as.data.frame(p) %>% mutate(biomarker = bio, model = "Single model")
  return(output)
}

#----------------------------------------------------------------------------------------------------
## Global model

get_pred2s <- function(out, bio, age, serotype, dat) {

  # Perform model
  if (out=="sev.only") {
    f <- paste(out, "~ (VCAM + SDC + Ang + IL8 + IP10 + IL1RA + CD163 + TREM + Fer + CRP) * Age * Serotype + uVCAM + uAng + uIL8 + uIP10 + uTREM")
  } else {
    f <- paste(out, "~ (rcs(VCAM, parms=kB(VCAM)) + rcs(SDC, parms=kB(SDC)) + rcs(Ang, parms=kB(Ang)) + rcs(IL8, parms=kB(IL8)) + rcs(IP10, parms=kB(IP10)) + rcs(IL1RA, parms=kB(IL1RA)) + rcs(CD163, parms=kB(CD163)) + rcs(TREM, parms=kB(TREM)) + rcs(Fer, parms=kB(Fer)) + rcs(CRP, parms=kB(CRP))) * rcs(Age, kB(Age)) +
                       (rcs(VCAM, parms=kB(VCAM)) + rcs(SDC, parms=kB(SDC)) + rcs(Ang, parms=kB(Ang)) + rcs(IL8, parms=kB(IL8)) + rcs(IP10, parms=kB(IP10)) + rcs(IL1RA, parms=kB(IL1RA)) + rcs(CD163, parms=kB(CD163)) + rcs(TREM, parms=kB(TREM)) + rcs(Fer, parms=kB(Fer)) + rcs(CRP, parms=kB(CRP))) * Serotype + 
                       uVCAM + uAng + uIL8 + uIP10 + uTREM")
  }
  
  if (out=="sev.or.inte") {
    m <- lrm(formula(f), data=dat)
  } else {
    m <- robcov(lrm(formula(f), data=dat, weights=ipwsd, x=T, y=T))
  }
  
  # Get prediction
  text <- paste("Predict(m, ", bio, "=dat$", bio, ", Age=", age, ", Serotype='", serotype, "', ref.zero=T)", sep="")
  p <- eval(parse(text = text))
  output <- as.data.frame(p) %>% mutate(biomarker = bio, model = "Global model")
  return(output)
}


#----------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------
# Function to get results from models with the interaction with serotype

#----------------------------------------------------------------------------------------------------
## Single model

get_est1s <- function(out, bio, ref1, ref2, est, age, serotype, dat) {

  # Perform model
  if (bio %in% c("VCAM", "Ang", "IL8", "IP10", "TREM", "Vir")) {
    if (out=="sev.only") {
      f <- paste(out, " ~ ", bio, " * Age * Serotype2 + u", bio, sep="")
    } else {
      f <- paste(out, " ~ rcs(", bio, ", kB(", bio, ")) * rcs(Age, kB(Age)) + rcs(", bio, ", kB(", bio, ")) * Serotype2 + u", bio, sep="")
    }
  } else {
    if (out=="sev.only") {
      f <- paste(out, " ~ ", bio, " * Age * Serotype2", sep="")
    } else {
      f <- paste(out, " ~ rcs(", bio, ", kB(", bio, ")) * rcs(Age, kB(Age)) + rcs(", bio, ", kB(", bio, ")) * Serotype2", sep="")
    }
  }
  
  if (out=="sev.or.inte") {
    m <- lrm(formula(f), data=dat)
  } else {
    m <- robcov(lrm(formula(f), data=dat, weights=ipwsd, x=T, y=T))
  }
  
  # Get p-values for (1) main effect, (2) interaction, and (3) non-linear effect
  if (est=="p") {output <- anova(m)[1,3]}
  else {
    if (est=="p int") {output <- anova(m)[2,3]}
    else {
      if (est=="p ns") {output <- anova(m)[3,3]}
      else {
        if (est=="p int age") {output <- anova(m)[10,3]}
        else {
          if (est=="p int serotype") {output <- anova(m)[16,3]}
          else {
            
            # Get OR and 95% CI
            text <- paste("summary(m, ", bio, "=c(", ref1, ",", ref2, "), ", "Age=", age, ", Serotype2=", serotype, ")", sep="")
            s <- eval(parse(text = text))
            if (est %in% c("OR", "HR")) {output <- s[2,4]}
            else {
              if (est=="loCI") {output <- s[2,6]}
              else {output <- s[2,7]}
            }
          }
        }
      }
    }
  }
  return(output)
}

#----------------------------------------------------------------------------------------------------
## Global model

get_est2s <- function(out, bio, ref1, ref2, est, age, serotype, dat) {

  # Perform model
  if (out=="sev.only") {
    f <- paste(out, "~ (VCAM + SDC + Ang + IL8 + IP10 + IL1RA + CD163 + TREM + Fer + CRP) * Age * Serotype2 + uVCAM + uAng + uIL8 + uIP10 + uTREM")
  } else {
    f <- paste(out, "~ (rcs(VCAM, parms=kB(VCAM)) + rcs(SDC, parms=kB(SDC)) + rcs(Ang, parms=kB(Ang)) + rcs(IL8, parms=kB(IL8)) + rcs(IP10, parms=kB(IP10)) + rcs(IL1RA, parms=kB(IL1RA)) + rcs(CD163, parms=kB(CD163)) + rcs(TREM, parms=kB(TREM)) + rcs(Fer, parms=kB(Fer)) + rcs(CRP, parms=kB(CRP))) * rcs(Age, kB(Age)) + 
                       (rcs(VCAM, parms=kB(VCAM)) + rcs(SDC, parms=kB(SDC)) + rcs(Ang, parms=kB(Ang)) + rcs(IL8, parms=kB(IL8)) + rcs(IP10, parms=kB(IP10)) + rcs(IL1RA, parms=kB(IL1RA)) + rcs(CD163, parms=kB(CD163)) + rcs(TREM, parms=kB(TREM)) + rcs(Fer, parms=kB(Fer)) + rcs(CRP, parms=kB(CRP))) * Serotype2 + 
                       uVCAM + uAng + uIL8 + uIP10 + uTREM")
  }
  
  if (out=="sev.or.inte") {
    m <- lrm(formula(f), data=dat)
  } else {
    m <- robcov(lrm(formula(f), data=dat, weights=ipwsd, x=T, y=T))
  }
  
  mp <- anova(m)
  
  # Get p-values for (1) main effect, (2) interaction, and (3) non-linear effect
  if (est=="p") {output <- mp[which(rownames(mp)==paste(bio, "(Factor+Higher Order Factors)", sep="  ")), 3]}
  else {
    if (est=="p int") {output <- mp[which(rownames(mp)==paste(bio, "(Factor+Higher Order Factors)", sep="  "))+1, 3]}
    else {
      if (est=="p ns") {output <- mp[which(rownames(mp)==paste(bio, "(Factor+Higher Order Factors)", sep="  "))+2, 3]}
      else {
        if (est=="p int age") {output <- mp[which(rownames(mp)==paste(bio, " * Age  (Factor+Higher Order Factors)", sep="")), 3]}
        else {
          if (est=="p int serotype") {output <- mp[which(rownames(mp)==paste(bio, " * Serotype2  (Factor+Higher Order Factors)", sep="")), 3]}
          else {
            
            # Get OR and 95% CI
            text <- paste("summary(m, ", bio, "=c(", ref1, ",", ref2, "), ", "Age=", age, ", Serotype2=", serotype, ")", sep="")
            s <- eval(parse(text = text))
            s0 <- s[which(rownames(s)==bio)+1,]
            if (est %in% c("OR", "HR")) {output <- s0[4]}
            else {
              if (est=="loCI") {output <- s0[6]}
              else {output <- s0[7]}
            }
          }
        }
      }
    }
  }
  return(output)
}


#----------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------