smooth_DimPlot <- function(sobj,
                           reduction,
                           group.by=NULL,
                           pt.size=.2,
                           label=FALSE,
                           repel=FALSE,
                           colors.use=c('red','blue'),
                           min.dens=1.e-5,
                           legend.size=8,
                           label.size=3) {
  
  require(MASS)
  require(akima)
  
  emb <- Embeddings(sobj, reduction)
  
  get_interpolated_dens <- function(emb, take) {
    
    dens <- kde2d(emb[take,1],emb[take,2])
    dens.df <- data.frame(with(dens,expand.grid(x,y)),as.vector(dens$z))
    colnames(dens.df) <- c('x','y','dens')
    with(dens.df, akima::interpp(x,y,dens,xo=emb[,1],yo=emb[,2],linear=F,extrap=T)$z)
  } 
  
  if (!is.null(group.by)) {
    gr <- droplevels(as.factor(FetchData(sobj, group.by)[,1]))
    levs <- levels(gr)
    if (length(levs) > 2) 
      warning("warning: more than 2 levels of grouping factors")
    
    cells.1 <- gr==levs[1]
    cells.2 <- gr==levs[2]
    
    dens.1 <- pmax(get_interpolated_dens(emb,cells.1),min.dens)
    dens.2 <- pmax(get_interpolated_dens(emb,cells.2),min.dens)
    
    log2ratio <- log2(dens.2/dens.1)
    
    if (all(levs %in% names(colors.use))) {
      col.low <- colors.use[levs[1]]
      col.high <- colors.use[levs[2]]
    } else {
      col.low <- colors.use[1]
      col.high <- colors.use[2]
    }
    
    df <- data.frame(x=emb[,1],
                     y=emb[,2],
                     log2ratio=log2ratio)
    
    pl <- ggplot(df,aes(x=x,y=y,colour=log2ratio)) +
      geom_point(size=pt.size, shape=16) + 
      scale_colour_gradient2(low=col.low,mid='gray',high=col.high) + 
      theme_bw()+
      theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
      coord_fixed(ratio=1)+
      guides(colour=guide_colorbar(title=paste0('log2ratio ',levs[2],' vs ',levs[1]),
                                   title.position='right',
                                   title.hjust=.5,
                                   title.theme=element_text(angle=90,size=legend.size),
                                   barwidth=.5))
  } else {
    dens <- get_interpolated_dens(emb, TRUE)
    df <- data.frame(x=emb[,1],
                     y=emb[,2],
                     dens=dens)
    
    pl <- ggplot(df,aes(x=x,y=y,colour=dens)) +
      geom_point(size=pt.size, shape=16) + 
      scale_colour_gradient(low=colors.use[1],high=colors.use[2]) +
      theme_bw()+
      theme(axis.ticks=element_blank(), panel.grid.minor=element_blank())+
      coord_fixed(ratio=1)+
      guides(colour=guide_colorbar(title=paste0('cell density'),
                                   title.position='right',
                                   title.hjust=.5,
                                   nbin=0,
                                   title.theme=element_text(angle=90,size=legend.size),
                                   barwidth=.5))
  }
  if (label) {
    df$idents <- Idents(sobj)
    mdf <- df %>% 
      dplyr::group_by(idents) %>%
      dplyr::summarise(x=median(x),
                       y=median(y))
    
    if (repel) {
      require(ggrepel)
      pl + geom_text_repel(data=mdf, aes(x=x,y=y,label=idents), colour='black', size=label.size)
    } else {
      pl + geom_text(data=mdf, aes(x=x,y=y,label=idents), colour='black', size=label.size)
    }
  } else {
    pl
  }
}
