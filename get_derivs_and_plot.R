### add derivs and plot function
font_size <- 18
theme_set(theme_classic(base_family = "sans",base_size = font_size))
theme_replace(axis.text=element_text(colour = "black",size = font_size))
line_size <- 1.5
point_size <- 2

### function to extract derivative, confidence interval, significance, and plot for GAMs ###
get_derivs_and_plot <- function(modobj,smooth_var,low_color=NULL,hi_color=NULL,xlabel=NULL,border_colors=NULL,flat_fill=FALSE){
  this_font_size = font_size
  if (is.null(low_color)){low_color = "white"}
  if (is.null(hi_color)){hi_color = "grey20"}

  derv<-derivatives(modobj,term=sprintf('s(%s)',smooth_var),partial_match = TRUE,level=0.99)
  derv<- derv %>%
    mutate(sig = !(0 >lower & 0 < upper))
  derv$sig_deriv = derv$derivative*derv$sig

  if ("unordered"%in%names(modobj$model)) {
    derv$fac_levels = factor(
      sapply(
        derv$smooth,
        FUN = function(x){str_match(x,levels(modobj$model$unordered))[!is.na(str_match(x,levels(modobj$model$unordered)))]}
      ),
      levels = levels(modobj$model$unordered)
    )
  } else {
    derv$fac_levels = factor(x=smooth_var,levels = smooth_var)
  }


  if (!is.null(border_colors)) {
    # get outline colors
    outline <- data.frame(
      fac_levels = factor(levels(derv$fac_levels),levels=levels(derv$fac_levels)),
      outline_color = border_colors
    )
  }

  if (flat_fill==TRUE) {
    # dplot <- ggplot(data=derv) + geom_tile(aes(x = data, y = .5, fill = sig))+
    #   scale_fill_manual(values = c("TRUE"="gray80","FALSE"="white"))+
    #   facet_grid(rows="fac_levels")+
    #   theme(panel.spacing = unit(-.01,"cm"))
    dplot<-ggplot(data=derv) + geom_tile(aes(x = data, y = .5, fill = abs(sig_deriv)))+
      scale_fill_gradient(low = "white",high = "black")+
      facet_grid(rows="fac_levels")+
      theme(panel.spacing = unit(-.01,"cm"))
  } else{
    dplot <- ggplot(data=derv) + geom_tile(aes(x = data, y = .5, fill = sig_deriv))+
      facet_grid(rows="fac_levels")+
      theme(panel.spacing = unit(-.01,"cm"))
    dplot <- dplot +
      scale_fill_gradient2(low = "#377eb8", midpoint = 0, mid = "white",
                           high = "#e41a1c",limits = c(min(derv$sig_deriv),max(derv$sig_deriv)))
  }

  # report numerical results
  factor_levels <-levels(derv$fac_levels)
  for (f in 1:length(factor_levels)) {
    this_level = factor_levels[f]
    this_derv <- derv %>%
      filter(str_detect(string = smooth,pattern = this_level))
    sig_ranges <- this_derv %>%
      mutate(blocks = rleidv(sig))%>%
      mutate(is_pos = derivative>0)%>%
      filter(sig==TRUE)%>%
      group_by(blocks)%>%
      summarize(min = min(data),max=max(data),positive = is_pos[1])

    cat(sprintf("\nSig change in %s:\n",this_level))
    if (length(sig_ranges$min)==0) {
      cat("No change\n")
    }
    else {
      for (r in 1:length(sig_ranges$blocks)) {
        if (sig_ranges$positive[r]==TRUE) {
          cat(sprintf("Increasing from %1.3f to %1.3f\n",sig_ranges$min[r],sig_ranges$max[r]))
        } else {
          cat(sprintf("Decreasing from %1.3f to %1.3f\n",sig_ranges$min[r],sig_ranges$max[r]))
        }
      }
    }

  }
  # }

  if (!is.null(xlabel)) {
    dplot <- dplot +
      labs(x = xlabel,fill = sprintf("\u0394%s",smooth_var))
  } else{
    dplot <- dplot +
      labs(x = smooth_var,fill = sprintf("\u0394%s",smooth_var))
  }
  dplot <- dplot +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.text.x = element_text(size = this_font_size),
          axis.line = element_blank(),
          axis.ticks.y = element_blank(),
          strip.background = element_blank(),
          strip.text.y = element_blank(),
          text = element_text(size=this_font_size),
          legend.text = element_text(size = this_font_size),
          axis.title = element_text(size = this_font_size),
          legend.key.width = unit(4,"cm"),
          legend.position = "bottom",
          legend.background = element_blank(),
          plot.margin = unit(c(0, 0, 0, 0), "cm"))+
    guides(fill = guide_colorbar(reverse = F,direction = "horizontal",title.position = "left"))
  if (!is.null(border_colors)) {
    dplot <- dplot +
      geom_rect(data=outline,
                aes(ymin=-0.001,ymax=1.001,xmin=min(derv$data)-max(derv$data)/1000,xmax=max(derv$data)+max(derv$data)/1000,color=border_colors),
                fill="white",alpha = 0,size=2)+
      scale_color_identity()

  } else{
    dplot <- dplot+geom_rect(aes(ymin=0,ymax=1,xmin=min(data),xmax=max(data)),color="black",fill="white",alpha = 0)
  }

  return(dplot)
}

