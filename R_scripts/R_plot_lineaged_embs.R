
vars <- c("AB_rel","Size","height","AR")
ylab <- c("Relative AB size","Embryo Size [um2]","Compression [um]", "Length / Width")

Embs$G1 <- Embs$Group
levls <- c("WT","Ctrl","Alive","Dead")
levels(Embs$G1) <- levls

Tuk <- function(var, groupvar, data) { #TukeyHSD returnig compact letter display in addition for given variable
  forml <- eval(substitute(v ~ gr, list(v = as.name(var), gr = as.name(groupvar))))
  mod <- TukeyHSD(aov(forml, data=data))
  pv <- mod[[groupvar]][,4]
  comp <- multcompLetters(pv)$Letters
  levls <- levels(data[,groupvar])
  comp <- comp[match(levls, names(comp))]
  output <- list(formula = forml, model = mod, comp = comp)
  attributes(output$formula) <- NULL
  return(output)
}

#for a given variable, which combinations are significant? -> feed them into a list that will be queried in the plotting function
ll <- list()
  for (v in vars) {
    comparisons <- Tuk(v, "G1", data=Embs) # perform multiple comparisons
    ymax <- as.vector(by(Embs[,v], Embs$G1, max, simplify = T)) #max for each group
    
    capture.output(comparisons, file = paste0(outDir,today,"_lineaged_Embs_phys_vars.txt"), append=T)
    
    #write.csv(pv, paste0(outDir,today,"_PhysicalParams_lineaged_Embs.csv"))
    ll[[v]] <- list(ymax=ymax, grcomp = comparisons$comp)
  }
  
pl <- ggboxplot(
          Embs[,c(vars,"G1")], "G1", vars,
          xlab=FALSE,
          ylab=ylab,
          title = c("AB size","Embryo Size","Compression", "Aspect ratio"),
          font.label = list(size = 8, color = "black"),
          add="jitter",
          add.params=list(size=1.5, shape=16, color="G1"),
          palette = pal,
          color="black",
          outlier.shape=NA,
          legend="none",
          bxp.errorbar=T,
          bxp.errorbar.width=0.2,
          error.plot = "errorbar"
          )
y.text=c(0.4,1100,16,1.35)

npl <- lapply(pl, function(x){
  a <- parent.frame()$i;
  x + theme(text = element_text(size=10),
          axis.text.y = element_text(size=8, angle=90, hjust = 0.5),
          axis.text.x = element_text(angle=90, hjust = 1, vjust = 0.5)
          )+
    annotate("text", x=1:4, y=ll[[a]]$ymax*1.05, label=ll[[a]]$grcomp) +
    annotate("segment", x=1.75,xend=4.25, y=y.text[a]-0.02*y.text[a], yend=y.text[a]-0.02*y.text[a]) +
    annotate("text",x=3,y=y.text[a], label="lin-5(ev571)", size=3, fontface = 'italic')
          })

ggarrange(plotlist=npl, ncol=4)
ggsave(paste0(outDir,today,"_PhysicalParameters_lineaged_Embs.pdf"), width =6 ,height=3)

rm(ll,pl, npl, vars, ylab, y.text, comparisons)
