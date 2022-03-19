library(tidyverse)
library(readxl)
library(resample)
library(ggsignif)
library(RColorBrewer)
library(ggpubr)
library(drawProteins)
library(ggcyto)

# functions necessary for confidence intervals -------
cibox<-function(x){
  b<-resample::bootstrap(x,R=30000,statistic=mean);
  m = b$observed;
  data.frame(lower=CI.percentile(b)[1],
             upper=CI.percentile(b)[2],
             ymin=CI.percentile(b)[1],
             ymax=CI.percentile(b)[2],
             middle=m)
}

permtest <- function(x,y,permutations=20000,statistic=mean) {
  b <- resample::permutationTest2(data=x, data2 = y,
                                  R=permutations,statistic=statistic)
  list( p.value = b$stats$PValue )
}
permmean <- function(x) {
  b = resample::bootstrap(x,R=5000,statistic=mean)$observed;
  data.frame(y = b,
             ymin = b,
             ymax = b)
}

# -------

setwd("~/Desktop/Humanized Mice Project (Jen)/")

engraftment<-read_csv("STAT1_EngraftmentTable1.csv") %>% 
  rename(Sample=...1) %>% 
  dplyr::filter(!grepl("SS_CD4",Sample),Sample!="Mean",Sample!="SD",Sample!="unstained001.fcs")

engraftment<-engraftment %>% separate(Sample,into=c("Sample",'Sample2','Sample3'),"_")

engraftment<-engraftment %>% mutate(Sample=case_when(
  grepl("Control",Sample) ~ paste0("C",Sample3),
  grepl("Exp",Sample) ~ paste0("E",Sample3),
  T ~ Sample
)) %>% select(-Sample2,-Sample3)
# Random quality checks.
unique(engraftment$Sample)
colnames(engraftment)

# Renaming
engraftment<-engraftment %>% rename(Monocyte_FreqTotal="FSC, SSC subset/Monocytes | Freq. of Total (%)",
                        CD3Tcells_FreqTotal="FSC, SSC subset/CD3_Tcells | Freq. of Total (%)" ,
                        CD19Bcells_FreqTotal="FSC, SSC subset/CD19_Bcells | Freq. of Total (%)" ,
                        CD4Tcells_FreqParent="FSC, SSC subset/CD3_Tcells/CD4_Tcells | Freq. of Parent (%)" ,
                        CD8Tcells_FreqParent="FSC, SSC subset/CD3_Tcells/CD8_Tcells | Freq. of Parent (%)" ,
) %>% select(-...9,-"FSC, SSC subset/CD3_Tcells/CD4_Tcells | Freq. of Total (%)",
             -"FSC, SSC subset/CD3_Tcells/CD8_Tcells | Freq. of Total (%)") %>% 
  dplyr::filter(Sample!="C10",Sample!="C11",Sample!="E18")

engraftment<-engraftment %>% 
  pivot_longer(-c(Sample)) %>% separate(name,into=c("Cell","Measure"),"_") %>% 
  #pivot_wider(values_from="value",names_from="name") %>% 
  mutate(Group=case_when(
    grepl("C",Sample) ~ "Control",
    grepl("E",Sample) ~ "Experimental"
  ))

engraftment<-engraftment %>% mutate(Cell=case_when(
  Cell=="Monocyte" ~ "Monocytes",
  Cell=="CD3Tcells" ~ "Total CD3 T Cells",
  Cell=="CD4Tcells" ~ "CD4 T Cells",
  Cell=="CD19Bcells" ~ "B Cells",
  Cell=="CD8Tcells" ~ "CD8 T Cells",
  T ~ Cell
))
engraftment$Cell<-factor(engraftment$Cell,levels=c("B Cells",
                                                   "Total CD3 T Cells",
                                                   "Monocytes",
                                                   "CD8 T Cells",
                                                   "CD4 T Cells"))

engraftment<-engraftment %>% mutate(Sample2=case_when(
  grepl("C",Sample) ~ "Control",
  grepl("17",Sample) ~ "E17",
  grepl("14",Sample) ~ "E14",
  T ~ Group
))

CD4CD8Ratio<-engraftment %>% 
  dplyr::filter(Cell=="CD4 T Cells"|Cell=="CD8 T Cells") %>% 
  pivot_wider(names_from=Cell,values_from=value) %>% 
  rename(CD4=`CD4 T Cells`,CD8=`CD8 T Cells`) %>% 
  mutate(CD4CD8Ratio=CD4/CD8) %>% 
  pivot_longer(!c(Sample,Sample2,Measure,Group)) %>% 
  dplyr::filter(name=="CD4CD8Ratio") %>% 
  rename(Cell=name)

engraftmentPlot<-rbind(engraftment,CD4CD8Ratio)

engraft<-ggplot(engraftmentPlot,aes(x=Group,y=value,color=Sample2,fill=Group))+
  geom_jitter(alpha=.9,position=position_jitterdodge())+
  theme_minimal()+
  scale_color_manual(values=c("#676664","#6EDED4","#6EDEF8","#FABA84"))+
  facet_wrap(.~Cell,scales="free_y")+
  ylab("% Frequency (All of Total, Except CD4/CD8 T Cells)")+
  #geom_bar(position = "dodge", stat = "summary", fun=mean,alpha=.1)+
  scale_fill_manual(values=c("#676664","#FABA84"))+
  stat_summary(alpha=1,fun.data=cibox, position=position_dodge(0.95),geom='errorbar',color="black",
               width=.2) +
  geom_signif(color="black",comparisons = list(c("Control","Experimental")),
              test=permtest, test.args = (permutations=20000),
              step_increase = 0.12, textsize=2.5,
              map_signif_level=T,
              margin_top=0)+
  theme(strip.text.x = element_text(size = 8),legend.position = "none",
        axis.text.x = element_text(size=7))+
  NULL

setwd("/Users/humzakhan/Desktop/Humanized\ Mice\ Project\ (Jen)/Figures")

ggsave("Figure3_Engraftment.png",plot=engraft,width=175,height=135,units="mm",scale=1)

#figure 2
setwd("/Users/humzakhan/Desktop/Humanized\ Mice\ Project\ (Jen)")
synthego<-read_xlsx("IndelPercent_Synthego.xlsx")
synthego[3,1]<-"Control"

ggplot(synthego,aes(x=Sample,y=IndelPercent))+
  geom_col(alpha=.45,fill="#73B5FE",color="grey")+
  ylab("Percent Indel/Knockout Score")+
  ggtitle("Indel Percentages in \nSTAT1 CRISPR Cells")+
  theme_bw()+
  ylim(0,100)+
  theme(text = element_text(size=7))

setwd("/Users/humzakhan/Desktop/Humanized\ Mice\ Project\ (Jen)/Figures")
ggsave("Figure2a_IndelPercentages.png",width=55,height=55,units="mm",scale=1,dpi=5000)

# gene annotation

stat1<-drawProteins::get_features("P42224")

stat1df<-drawProteins::feature_to_dataframe(stat1)

p <- draw_canvas(stat1df)
p <- draw_chains(p, stat1df,fill="lightsteelblue1",outline = "grey",label_chains=F)
p <- draw_domains(p, stat1df)
p <- draw_regions(p, stat1df)
p <- draw_motif(p, stat1df)
p <- draw_phospho(p, stat1df, size = 8) 

p <- p + theme_bw(base_size = 5) + theme_bw() + theme(axis.text.y = element_blank(),
                                                       axis.ticks.y=element_blank())+
  geom_vline(xintercept = 268,linetype="dashed")+theme(legend.position = "none")+
  theme(text=element_text(size=5))

p

ggsave("Figure2b_STAT1_Protein.png",width=110,height=55,units="mm",scale=1,dpi=2500)

# figure 4 ------

setwd("/Users/humzakhan/Desktop/Humanized\ Mice\ Project\ (Jen)/2022_STAT1_CSVs/")

cells<-as_tibble(NULL)
for (i in list.files()){
  sample<-str_split(i,"_")[[1]][2]
  stim<-str_split(i,"_")[[1]][3]
  cell<-str_split(i,"_")[[1]][4]
  print(i)
  exprs<-read_csv(i)
  # doesnt work for SS (as it shouldnt)
  colnames(exprs)[grepl("PE",colnames(exprs))&!grepl("PE-Cy7",colnames(exprs))]<-"CD4"
  colnames(exprs)[grepl("AF647",colnames(exprs))]<-"TotalSTAT1"
  colnames(exprs)[grepl("AF488",colnames(exprs))]<-"pSTAT1"
  colnames(exprs)[grepl("APC-Cy7",colnames(exprs))|grepl("PE-Cy7",colnames(exprs))]<-"CD19"
  colnames(exprs)[grepl("BV421",colnames(exprs))]<-"CD14"
  
  exprs$Sample<-sample
  exprs$Stim<-stim
  exprs$Cell<-cell
  exprs <- exprs %>% select(pSTAT1,CD4,Sample,Stim,Cell,TotalSTAT1,CD14,CD19)
  cells<-rbind(exprs,cells)
}


cells<-cells %>% mutate(Cell=case_when(
  Cell=="Q3 CD19 APC-Cy7+ , CD4 PE-.csv" ~ "B Cell", 
  Cell=="Q1 CD19 APC-Cy7- , CD4 PE+.csv"~ "CD4 T Cell",
  Cell=="Monocytes.csv" ~ "Monocyte"
))

cells<-cells %>% mutate(Stim=case_when(
  grepl("IgG",Stim) ~ "IgG Isotype", 
  grepl("alpha",Stim) ~ "IFNa", 
  grepl("gamma",Stim) ~ "IFNg", 
  grepl("unstim",Stim) ~ "Unstimulated", 
))

cells<-cells %>% mutate(Group=case_when(
  grepl("E",Sample) ~ "Experimental",
  grepl("C",Sample) ~ "Control"
))

cells<-cells %>% dplyr::filter(Sample!="C10",Sample!="C11",Sample!="E18")

cells<-cells %>% group_by(Stim,Cell) %>% mutate(pSTAT1Norm=pSTAT1/median(pSTAT1),
                                                pSTAT1NormtSTAT=pSTAT1/median(TotalSTAT1),
                                             TotalSTAT1Norm=TotalSTAT1/median(TotalSTAT1))

adjFact<-(1+abs(min(cells$pSTAT1Norm)))
adjFactptSTAT<-(1+abs(min(cells$pSTAT1NormtSTAT)))


cells<-cells %>% mutate(pSTAT1Norm=pSTAT1Norm+adjFact,
                        pSTAT1NormtSTAT=pSTAT1NormtSTAT+adjFactptSTAT)

means<-cells %>%  group_by(Stim,Cell,Group) %>% summarise(tSTAT1Med=median(TotalSTAT1),
                                                          TotalSTAT1NormMed=median(TotalSTAT1Norm),
                                                         pSTAT1Med=median(pSTAT1),
                                                         pSTAT1NormMed=median(pSTAT1Norm),
                                                         pSTAT1NormtSTATMed=median(pSTAT1NormtSTAT)
                                                         )

# TSTAT
cells<-cells %>% 
  mutate(Condition=case_when(
    grepl("C",Sample) ~ "Control",
    T ~  Sample
  ))

ggplot(cells %>% dplyr::filter(grepl("14|17|C",Sample),
                               Stim=="Unstimulated"),
       aes(x=TotalSTAT1,fill=Condition))+
  facet_grid(.~Cell)+
  geom_density(alpha=.5)+
  theme_bw()+
  xlab("Total STAT1 Protein")+
  ylab("Density")+
  theme(strip.text.x = element_text(size=7),
        strip.background = element_rect("#F5F5F5"),
        legend.key.height= unit(4, 'mm'),
        legend.key.width= unit(4, 'mm'))+
  #scale_x_flowjo_biexp()+
  scale_x_log10()+
  geom_vline(data=means %>% dplyr::filter(Stim=="Unstimulated",grepl("Control",Group)),
             aes(xintercept=tSTAT1Med),
             linetype="dashed")+
  annotation_logticks(sides = "b",size=.1)+  
  scale_fill_manual(values=c("#FFE9B3","#6EDED4","#6EDEF8"))+
  NULL

# regular looks just as good as median normalized
setwd("/Users/humzakhan/Desktop/Humanized\ Mice\ Project\ (Jen)/Figures")
ggsave("Figure4a_TotalSTAT1_Protein.png",width=175,height=55,units="mm",scale=1,dpi=2500)

# ggplot(cells %>% dplyr::filter(!grepl("14|17",Sample),
#                                Stim=="Unstimulated"),
#        aes(x=TotalSTAT1,fill=PlotGroup))+
#   facet_wrap(.~Cell,scales="free_y",ncol=3)+
#   geom_density(alpha=.5)+
#   theme_bw()+
#   #scale_x_log10()+
#   scale_x_flowjo_biexp()+
#   xlab("Total STAT1 Protein")+
#   ylab("Density")+
#   theme(strip.text.x = element_text(size=15))+
#   geom_vline(data=means %>% dplyr::filter(Stim=="Unstimulated",grepl("Control",Group)),
#              aes(xintercept=tSTAT1Med),
#              linetype="dashed")+
#   scale_fill_brewer(palette="BuPu")+
#   NULL

# pSTAT

ggplot(cells %>% dplyr::filter(grepl("14|17|C",Sample),pSTAT1Norm<10),
       aes(x=pSTAT1Norm,fill=Condition))+
  facet_wrap(Stim~Cell,scales="free_y",ncol=3)+
  geom_density(alpha=.4)+
  theme_bw()+
  scale_x_log10()+
  #scale_x_flowjo_biexp()+
  xlab("pSTAT1 (Normalized to Median)")+
  ylab("Density")+
  theme(strip.text.x = element_text(size=10),strip.background = element_rect("#F5F5F5"))+
  geom_vline(data=means %>% dplyr::filter(grepl("Control",Group)),
             aes(xintercept=pSTAT1NormMed),
             linetype="dashed")+
  scale_fill_manual(values=c("#FFE9B3","#6EDED4","#6EDEF8"))+
  NULL

ggsave("Figure4b_pSTAT1_StimAssays_NORM.png",width=175,height=150,units="mm",scale=1,dpi=3000)

 # ggplot(cells %>% dplyr::filter(!grepl("14|17",Sample)),
 #       aes(x=pSTAT1Norm,fill=Condition))+
 #  facet_wrap(Stim~Cell,scales="free_y",ncol=3)+
 #  geom_density(alpha=.35)+
 #  theme_bw()+
 #  scale_x_log10()+
 #  #scale_x_flowjo_biexp()+
 #  xlab("pSTAT1 (Normalized to Median)")+
 #  ylab("Density")+
 #  theme(strip.text.x = element_text(size=15))+
 #  geom_vline(data=means %>% dplyr::filter(grepl("Control",Group)),
 #             aes(xintercept=pSTAT1NormMed),
 #             linetype="dashed")+
 #  scale_fill_brewer(palette="BuPu")+
 #  NULL
rel3<-cells %>% dplyr::filter(grepl("14|17|C",Sample),
                              !grepl("IgG",Stim),
                              pSTAT1NormtSTAT<10)

rel3$Stim<-factor(rel3$Stim,levels = c("Unstimulated", "IFNa", "IFNg"))

ggplot(rel3,
       aes(x=pSTAT1NormtSTAT,fill=Condition))+
  facet_wrap(Stim~Cell,scales="free_y",ncol=3)+
  geom_density(alpha=.6)+
  theme_bw()+
  scale_x_log10()+
# scale_x_flowjo_biexp()+
  xlab("pSTAT1NormtSTAT")+
  ylab("Density")+
  theme(strip.text.x = element_text(size=10),strip.background = element_rect("#F5F5F5"))+
  geom_vline(data=means %>% dplyr::filter(grepl("Control",Group),!grepl("IgG",Stim)),
             aes(xintercept=pSTAT1NormtSTATMed),
             linetype="dashed")+
  scale_fill_manual(values=c("#FFE9B3","#6EDED4","#6EDEF8"))+
  NULL

ggsave("Figure4b_pSTAT1_StimAssays_NORMtotSTAT.png",width=175,height=150,units="mm",scale=1,dpi=3000)

rel4<-cells %>% dplyr::filter(grepl("14|17|C",Sample))

rel4$Stim<-factor(rel4$Stim,levels=c("IgG Isotype","Unstimulated","IFNa","IFNg")) %>% f

monoMeans<-cells %>% dplyr::filter(grepl("14|17|C",Sample),Cell=="Monocyte") %>% 
  group_by(Stim,Cell,Condition) %>% 
  summarise(pSTAT1Mean=mean(pSTAT1))

ggplot(cells %>% dplyr::filter(grepl("14|17|C",Sample),Cell=="Monocyte"),
       aes(x=pSTAT1,fill=Stim))+
  facet_wrap(Condition~.,scales="free_y",nrow=3)+
  geom_density(alpha=.4)+
  theme_bw()+
  #scale_x_log10()+
  scale_x_flowjo_biexp()+
  xlab("pSTAT1")+
  ylab("Density")+
  theme(strip.text.x = element_text(size=10),strip.background = element_rect("#F5F5F5"))+
  geom_vline(data=monoMeans,
             aes(xintercept=pSTAT1Mean),
             linetype="dashed")+
 scale_fill_manual(values=c("#B452F7","#526FF7","#FAD7D7","#F4FFB5"))+
  NULL

cd4Means<-cells %>% dplyr::filter(grepl("14|17|C",Sample),Cell=="CD4 T Cell") %>% 
  group_by(Stim,Cell,Condition) %>% 
  summarise(pSTAT1Mean=mean(pSTAT1))

ggplot(cells %>% dplyr::filter(grepl("14|17|C",Sample),Cell=="CD4 T Cell"),
       aes(x=pSTAT1,fill=Stim))+
  facet_wrap(Condition~.,scales="free_y",nrow=3)+
  geom_density(alpha=.4)+
  theme_bw()+
  #scale_x_log10()+
  scale_x_flowjo_biexp()+
  xlab("pSTAT1")+
  ylab("Density")+
  theme(strip.text.x = element_text(size=10),strip.background = element_rect("#F5F5F5"))+
  geom_vline(data=cd4Means,
             aes(xintercept=pSTAT1Mean),
             linetype="dashed")+
  scale_fill_manual(values=c("#B452F7","#526FF7","#FAD7D7","#F4FFB5"))+
  NULL

#ggsave("Figure4b_pSTAT1_StimAssays_NORM.png",width=175,height=150,units="mm",scale=1,dpi=3000)

# 
# ggplot(cells %>% dplyr::filter(!grepl("14|17|C",Sample)),
#        aes(x=pSTAT1,fill=PlotGroup))+
#   facet_wrap(Stim~Cell,scales="free_y",ncol=3)+
#   geom_density(alpha=.5)+
#   theme_bw()+
#   xlab("pSTAT1")+
#   ylab("Density")+
#   theme(strip.text.x = element_text(size=10),strip.background = element_rect("#F5F5F5"))+
#   scale_x_flowjo_biexp()+
#   geom_vline(data=means %>% dplyr::filter(grepl("Control",Group)),
#              aes(xintercept=pSTAT1Med),
#              linetype="dashed")+
#   scale_fill_brewer(palette="BuPu")+
#   NULL


# pstat by tstat 
relevant<-cells %>% dplyr::filter(grepl("14|17|C",Sample),!grepl("IgG",Stim))
relevant$Stim<-factor(relevant$Stim,levels=c("Unstimulated","IFNa","IFNg"))
means$Stim<-factor(means$Stim,levels=c("Unstimulated","IFNa","IFNg"))

ggplot(relevant,aes(x=TotalSTAT1,y=pSTAT1,color=Condition))+
  facet_wrap(Cell~Stim,ncol=3)+
  scale_x_flowjo_biexp()+
  scale_y_flowjo_biexp()+
  geom_point(alpha=.3)+
  theme_bw()+
  xlab("Total STAT1")+
  ylab("pSTAT1")+
  theme(strip.text.x = element_text(size=10),strip.background = element_rect("#F5F5F5"))+
  scale_color_manual(values=c("#FFE9B3","#FBA8A8","#6EDEF8"))+
  geom_vline(data=means %>% dplyr::filter(grepl("Control",Group),!grepl("IgG",Stim)),
             aes(xintercept=tSTAT1Med),
             linetype="dashed")+
  geom_hline(data=means %>% dplyr::filter(grepl("Control",Group),!grepl("IgG",Stim)),
             aes(yintercept=pSTAT1Med),
             linetype="dashed")+
  NULL

ggplot(cells %>% dplyr::filter(grepl("14|17|C",Sample)),
       aes(x=TotalSTAT1,y=pSTAT1,color=Condition))+
  #facet_wrap(Cell~Stim,ncol=4)+
  stat_density2d()+
  theme_bw()+
  #xlab("Total STAT1")+
 # ylab("pSTAT1")+
 # theme(strip.text.x = element_text(size=10),strip.background = element_rect("#F5F5F5"))+
  # scale_x_flowjo_biexp()+
  # scale_y_flowjo_biexp()+
  # scale_color_manual(values=c("#FFE9B3","#FBA8A8","#6EDEF8"))+
  # geom_vline(data=means %>% dplyr::filter(grepl("Control",Group)),
  #            aes(xintercept=tSTAT1Med),
  #            linetype="dashed")+
  # geom_hline(data=means %>% dplyr::filter(grepl("Control",Group)),
  #            aes(yintercept=pSTAT1Med),
  #            linetype="dashed")+
  NULL

cellmeds<-cells %>% 
  group_by(Sample,Stim,Cell,Condition) %>% 
  summarise_if(is.numeric,median)

cellmedsList<-cellmeds %>% split(list(cellmeds$Sample,cellmeds$Cell))

calculateFoldChange<-function(x) {
  pSTAT1US<-x %>% dplyr::filter(Stim=="Unstimulated") %>% pull(pSTAT1)  
  pSTAT1NormUS<-x %>% dplyr::filter(Stim=="Unstimulated") %>% pull(pSTAT1NormtSTAT)  
  
  IFNaFCNorm<-(x %>% dplyr::filter(Stim=="IFNa") %>% pull(pSTAT1NormtSTAT))/pSTAT1NormUS
  IFNaFC<-(x %>% dplyr::filter(Stim=="IFNa") %>% pull(pSTAT1))/pSTAT1US
  
  IFNgFCNorm<-(x %>% dplyr::filter(Stim=="IFNg") %>% pull(pSTAT1NormtSTAT))/pSTAT1NormUS
  IFNgFC<-(x %>% dplyr::filter(Stim=="IFNg") %>% pull(pSTAT1))/pSTAT1US
  
  return(tibble(Sample=x$Sample,Cell=x$Cell,IFNa_Norm=IFNaFCNorm,
                IFNa_NonNorm=IFNaFC,IFNg_Norm=IFNgFCNorm,
                IFNg_NonNorm=IFNgFC))
  
  
}

FCs<-map_dfr(cellmedsList,calculateFoldChange) %>% unique()

FCs<-FCs %>% pivot_longer(!c(Sample,Cell)) %>% 
  separate(name,into=c("Stim","Norm")) %>% mutate(Group=case_when(
    grepl("C",Sample)~"Control",
    grepl("E",Sample)~"Experimental"
  ))


ggplot(FCs,aes(x=Group,y=value,color=Group))+
  facet_wrap(Cell~Stim,ncol=2,scales="free_y")+
  geom_jitter(alpha=.3)+
  theme_bw()+
  stat_summary(alpha=1,fun.data=cibox, position=position_dodge(0.95),geom='errorbar',color="black",
               width=.2) +
  geom_signif(color="black",comparisons = list(c("Control","Experimental")),
              test=permtest, test.args = (permutations=20000),
              step_increase = 0.12, textsize=2.5,
              map_signif_level=T,
              margin_top=0)
  

ggsave("Figure_pSTAT_FoldChanges.png",width=175,height=150,units="mm",scale=1,dpi=3000)

  
  
  

