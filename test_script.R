---
  title: "Figure 1"
author: "Kverneland_AH"
date: "`r Sys.Date()`"
output: html_document
---
  ```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
library(tidyverse)
library(ggrepel)
library(readxl)
library(RColorBrewer)
library(corrplot)
library(ggVennDiagram)
library(ggpmisc)
library(ggpubr)


enrichment_color <- set_names( c("#FBB4AE", "#DECBE4", "#CCEBC5", "#B3CDE3", "#FED9A6"), c("Neat", "Top14", "MagNET", "EV20k", "AcidP"))


kba_color <- set_names(brewer.pal(10,"Set3"),c("ALB","SERPINA1","CRP","CST3", "HP","LPA","B2M","MB","TF","SHBG"))

brewer.pal(5,"Pastel1")
```

## Load data

```{r load data, echo=FALSE}
data_prp_full <- read_csv("outputs/fig1_prot_norm.csv", show_col_types = F) %>% 
  mutate(enrichment=factor(enrichment, levels=c("Neat", "Top14", "MagNET", "EV20k"))) 
data_prp_full %>% distinct(enrichment, material, gradient)

data_prp <- data_prp_full %>% 
  filter(enrichment!="Neat_Rep") %>% 
  filter(gradient=="60SPD")

annotation_samples <- read_csv("outputs/annotation_samples.csv", show_col_types = F)
annotation_prot <- read_csv("data/annotations_hpa_2024-09-11.csv", show_col_types = F) %>% 
  mutate(subcellular = case_when(deeploc_single=="Extracellular"~"Extracellular",
                                 deeploc_single=="Cell membrane"~"Cell membrane",
                                 deeploc_single!="Extracellular" | deeploc_single=="Cell membrane" ~ "Intracellular")
         
  )
kba_list <- read_csv("kba_list.csv", show_col_types = F)
kba_list$PG.Genes

annotation_prot %>% distinct(type)



```

## Counts
Identified protein counts

```{r count, echo=T}

data_prp_count <- data_prp %>%
  group_by(enrichment, gradient, samples, name_short) %>% 
  count()

prot_counts <- data_prp_count %>% 
  ggplot(aes(enrichment, n))+
  geom_col(aes(fill=enrichment), position = "dodge2", color="black")+
  scale_fill_manual(values=enrichment_color)+
  stat_summary(fun = "median", colour = "black", size = 4,
               geom = "text", aes(label = round(after_stat(y),0)),vjust=-.5)+
  theme_classic()+
  theme(legend.position = "none", axis.title.x = element_blank())+
  labs(title="Identified Protein Groups", x="Samples")+
  ylim(0,4700)
prot_counts  
ggsave("fig1/prot_counts.png",height=2.5, width=4)

```

## VENNN

```{r}
library(eulerr)
library(venn)
venn_list <- list(
  "Neat" = data_prp %>% filter(enrichment=="Neat") %>% distinct(PG.ProteinGroups) %>% pull(),
  "Top14" = data_prp %>% filter(enrichment=="Top14") %>% distinct(PG.ProteinGroups) %>% pull(),
  "MagNET" = data_prp %>% filter(enrichment=="MagNET") %>% distinct(PG.ProteinGroups)%>% pull(),
  "EV20k" = data_prp %>% filter(enrichment=="EV20k") %>% distinct(PG.ProteinGroups)%>% pull()
  
)
euler(venn_list)




library(VennDiagram)

venn<-VennDiagram::venn.diagram(venn_list,
                                filename = "fig1/venn_alt.png",
                                cex=2,
                                fontfamily="sans",
                                cat.cex=2,
                                cat.fontfamily = "sans",
                                fill=enrichment_color[1:4])


venn


```

## Subcellular

```{r subcellular}
data_prp_summary_n <- data_prp %>% 
  left_join(annotation_prot, "PG.ProteinGroups") %>% 
  count(enrichment, subcellular) %>% 
  group_by(enrichment) %>% 
  mutate(proportion=n/sum(n)*100) %>% 
  ggplot(aes(x=enrichment, y=proportion))+
  geom_col(aes(fill=subcellular), color="black")+
  theme_classic()+
  ylim(0,105)+
  labs(title="Subcellular location", y="Proportion of Protein Groups (%)", fill="Method")
data_prp_summary_n
ggsave("fig1/subcellular_n.png",height=2.5, width=4)

data_prp_subcell_int <- data_prp %>% 
  mutate(PG.ProteinGroups=sub("\\;.*", "", PG.ProteinGroups)) %>% 
  left_join(annotation_prot, "PG.ProteinGroups") %>% 
  filter(!is.na(annotation_alb)) %>% 
  mutate(annotation_func=if_else(annotation_alb %in% c("Albumin", "Immunoglobulins", "Coagulation"), annotation_alb, "Other")) %>% 
  group_by(enrichment, annotation_func) %>%
  summarize(LFQ_sum = sum(PG.Quantity, na.rm=T), .groups = "drop_last") %>%
  mutate(LFQ_prop = LFQ_sum/sum(LFQ_sum)*100) 

data_prp_subcell_int %>% 
  ggplot(aes(x=LFQ_prop, y=enrichment))+
  geom_col(aes(fill=enrichment, alpha=annotation_func),color="black")+
  scale_fill_manual(values=enrichment_color)+
  labs(title="Subcellular location (intensity)", y="Proportion of summed LFQ (%)", fill="Method", alpha="Subcellular")+
  theme_classic()+
  theme(legend.key.size = unit(0.4, 'cm'), axis.title.x = element_blank())
data_prp_subcell_int
ggsave("fig1/subcellular_lfq.png",height=2.5, width=4)
```


## Protein type

```{r type}
variable_test<- secretome_location

data_prp_protein_type <- data_prp %>% 
  left_join(annotation_prot, "PG.ProteinGroups") %>% 
  filter(!is.na(annotation_alb)) %>% 
  group_by(enrichment, annotation_alb) %>% 
  summarize(LFQ_sum=sum(PG.Quantity,na.rm=T), .group="drop_last") %>% 
  mutate(LFQ_prop= LFQ_sum/sum(LFQ_sum)*100) %>% 
  ggplot(aes(y=fct_reorder(annotation_alb, LFQ_prop),x=LFQ_prop))+
  geom_col(aes(fill=enrichment), position = "dodge2", color="black")+
  scale_fill_manual(values=enrichment_color)+
  theme_classic()+
  theme(axis.title.y = element_blank())+
  labs(title="Protein function", x="Proportion of intensity")
data_prp_protein_type
ggsave("fig1/prot_function.png", height=2.5,width=4)

data_prp_protein_type_n <- data_prp %>% 
  left_join(annotation_prot, "PG.ProteinGroups") %>% 
  filter(!is.na(annotation_alb)) %>% 
  group_by(enrichment, annotation_alb) %>% 
  count() %>% 
  ggplot(aes(y=fct_reorder(annotation_alb, n),x=n))+
  geom_col(aes(fill=enrichment), position = "dodge2")+
  scale_fill_manual(values=enrichment_color)+
  theme_classic()+
  theme(axis.title.y = element_blank())
data_prp_protein_type_n
```



## FDA biomarker target
```{r}
data_prp_summary_rank <- data_prp %>% 
  filter(!is.na(PG.Quantity)) %>% 
  group_by(enrichment, PG.ProteinGroups) %>% 
  summarise(n=n(),
            median_LFQ = median(PG.Quantity,na.rm=T),
            median_LFQ_log2=median(log2(PG.Quantity),na.rm=T),
            median_LFQ_log10 =median(log10(PG.Quantity),na.rm=T),
            cv=sd(PG.Quantity,na.rm=T)/mean(PG.Quantity,na.rm=TRUE)*100) %>% 
  arrange(enrichment, desc(median_LFQ)) %>% 
  mutate(rank_lfq=row_number()) %>% 
  arrange(cv) %>% 
  mutate(rank_cv=row_number()) %>%
  left_join(annotation_prot, "PG.ProteinGroups")


data_prp_summary_fda <- data_prp_summary_rank %>% 
  mutate(FDA_drug_target_binary = case_when(FDA_drug_target=="FDA_target"~"Yes", 
                                            is.na(FDA_drug_target) | FDA_drug_target!="FDA_target" ~"No")) %>% 
  mutate(FDA_drug_target_binary = factor(FDA_drug_target_binary, levels=c("Yes", "No")))


data_prp_summary_fda_neat<-data_prp_summary_fda %>% filter(enrichment=="Neat")

data_prp_summary_fda %>% 
  arrange(enrichment,log10(median_LFQ)) %>% 
  ggplot(aes(x=rank_lfq, y=median_LFQ_log10)) +
  geom_point(data=data_prp_summary_fda %>% filter(FDA_drug_target_binary=="No"),  aes(color=enrichment))+
  geom_point(data=data_prp_summary_fda %>% filter(FDA_drug_target_binary=="Yes"), aes(shape=FDA_drug_target_binary), alpha=0.5, color="darkgrey") +
  geom_text(data=data_prp_summary_fda %>% count(FDA_drug_target_binary,enrichment) %>% filter(FDA_drug_target_binary=="Yes"),
            aes(x=4000,y=5.5,label=paste0("n=",n)))+
  facet_wrap(~enrichment,ncol=1)+
  scale_shape_manual(values=c(22))+
  scale_color_manual(values=enrichment_color)+
  theme_classic()+
  labs(title="FDA drug targets", x="Rank in dataset", y="LFQ (log10)", color="Method", shape="FDA biomarker")
ggsave("fig1/fda_rank.png",height=3.5, width=5)
```


## Correlation vs neat

```{r}
data_prp_summary_neat <- data_prp %>% 
  filter(!is.na(PG.Quantity)) %>% 
  group_by(enrichment, PG.ProteinGroups) %>% 
  summarise(median_LFQ_log2=median(log2(PG.Quantity),na.rm=T),
            median_LFQ_log10 =median(log10(PG.Quantity),na.rm=T)) %>% 
  pivot_wider(id_cols=PG.ProteinGroups, names_from = enrichment, values_from = median_LFQ_log10) %>% 
  pivot_longer(3:5, names_to = "Method", values_to = "LFQ") %>% 
  left_join(annotation_prot, "PG.ProteinGroups")

data_prp_summary_neat %>%
  group_by(Method) %>% 
  filter(!is.na(Neat)) %>% 
  count(!is.na(Neat), !is.na(LFQ))




data_prp_summary_neat %>% 
  ggplot(aes(y=Neat, x=LFQ))+
  geom_point(aes(color=subcellular),size=1, alpha=0.7)+
  labs(title="Correlation to Neat plasma", y="LFQ in Neat (log10)", x="LFQ in method (log10")+
  stat_poly_eq(parse = TRUE)+
  geom_smooth(method="glm")+
  geom_text(data=data_prp_summary_neat %>% filter(!is.na(Neat) & !is.na(LFQ)) %>% count(Method), aes(label=paste0("n = ",n), x=7,y=17))+
  theme_classic()+
  facet_wrap(~Method)
ggsave("fig1/correlation_neat.png", height=2.5, width=7)



```




## KBA

```{r kba}


port_norm_kba_median <- port_prp_kba %>%  
  group_by(enrichment, material, PG.Genes) %>% 
  summarize(LFQ_median=median(log10(2^LFQ_log2)),
            LFQ_median_vsn=median(LFQ_vsn)) %>% 
  left_join(kba_list, "PG.Genes") %>% 
  left_join(annotation_prot, "PG.Genes") %>% 
  arrange(qc_rank) %>% 
  mutate(PG.Genes=factor(PG.Genes,levels=reorder(PG.Genes, desc(qc_rank)))) %>% 
  mutate(Below_LOD=if_else(PG.Genes %in% c("CRP", "LPA") , "Yes", "No")) %>% 
  mutate(Name_LOD=if_else(PG.Genes %in% c("CRP", "LPA") , paste0(PG.Genes,"*"), PG.Genes))

port_norm_kba_median %>% 
  pivot_wider(id_cols=c(PG.Genes, deeploc_single, qc_rank), names_from = c(material, enrichment), values_from="LFQ_median") %>% 
  arrange(qc_rank) 


port_norm_kba_median_LOD <- port_norm_kba_median %>%  filter(!PG.Genes %in% c("CRP", "LPA"))

port_norm_kba_median %>% 
  ggplot(aes(y=qc_conc, x=LFQ_median))+
  geom_smooth(data=port_norm_kba_median_LOD, method="lm", color="black")+
  stat_poly_eq(data=port_norm_kba_median_LOD, parse = TRUE)+
  geom_point(aes(color=Name_LOD))+
  geom_text_repel(data=port_norm_kba_median, aes(label=Name_LOD), size=3)+
  scale_y_log10(n.breaks=4, labels=c("ng/L","µg/L","mg/L","g/L"))+
  facet_grid(material~enrichment)+
  theme_classic()+
  theme(legend.key.height = unit(0.5, 'cm'))+
  labs(title="KBA correlation w/o CRP/LPA")
ggsave("fig1/KBA_correlation.png", height=2.5, width=7)
?stat_poly_eq
```




## TOP14

```{r top14}
top14 <- read_csv( "top14.csv")


```



# Supplementary

## Gradient
```{r KBA stability}

port_prp_kba <- data_prp %>% 
  filter(PG.Genes %in% kba_list$PG.Genes) %>% 
  mutate(Gene_name=factor(PG.Genes, 
                          levels=c("ALB", "SERPINA1","HP", "TF", "LPA", "SHBG", "B2M", "CST3", "CRP", "MB")))

kba_proteins_plot <-  port_prp_kba %>% 
  ggplot(aes(x=samples, y=LFQ_vsn))+
  geom_point(aes(color=Gene_name))+
  geom_line(aes(group=Gene_name, color=Gene_name))+
  facet_grid(material~enrichment)+
  theme_classic()+
  ylim(0,25)+
  theme(legend.key.height = unit(0.5, 'cm'))+
  scale_x_continuous(n.breaks=8)+
  labs(title="Protein quantity across replicates", x="Workflow replicates", y="LFQ (log2)", color="Protein name")
kba_proteins_plot
ggsave("fig1/kba_protein_stability.png", height=2.5, width=5)


```




```{r sup gradient saxvol}
prot_norm_all <- read_csv("outputs/fig2_prot_norm_all.csv")


prp_full_count <- prot_norm_all %>% 
  filter(enrichment!="Neat_Rep") %>% 
  group_by(enrichment, gradient, samples, sax_vol, name_short) %>% 
  count()

prp_full_count_plot <- prp_full_count %>% 
  filter(sax_vol!="SAX40" | is.na(sax_vol)) %>% 
  ggplot(aes(enrichment, n))+
  geom_col(aes(fill=enrichment), position = "dodge2", color="black")+
  scale_fill_manual(values=enrichment_color)+
  stat_summary(fun = "median", colour = "black", size = 4,
               geom = "text", aes(label = round(after_stat(y),0)),vjust=-0.5)+
  theme_classic()+
  theme(legend.position = "none")+
  labs(title="Gradient comparison", x="Samples")+
  facet_wrap(~gradient)+
  ylim(0,4900)
prp_full_count_plot  
ggsave("fig1/sup_gradient_counts.png",height=2.5, width=5)


```


## Complement
```{r type}
complement <- data_prp %>% 
  left_join(annotation_prot %>% dplyr::select (-PG.Genes), "PG.ProteinGroups") %>% 
  filter(type=="COMPLEMENT" | secretome_function=="Complement pathway") %>% 
  distinct(PG.Genes, PG.ProteinGroups, PG.ProteinDescriptions)
write_csv(complement, "outputs/complement.csv")

data_prp_complement <- data_prp %>% 
  group_by(enrichment,PG.ProteinGroups) %>% 
  summarize(median_LFQ = median(LFQ_log2,na.rm=T)) %>%
  mutate(rank=-rank(median_LFQ)) %>% 
  left_join(annotation_prot, "PG.ProteinGroups")

data_prp_complement %>%   
  ggplot(aes(x=rank, y=median_LFQ))+
  geom_point(aes(fill=enrichment), shape=21, color="lightgrey")+
  geom_point(data=data_prp_complement %>%   filter(type=="COMPLEMENT" | secretome_function=="Complement pathway"),
             aes(fill=enrichment), shape=21, color="black")+
  geom_text_repel(data=data_prp_complement %>%   filter(type=="COMPLEMENT" | secretome_function=="Complement pathway"),
                  aes(label=PG.Genes), size=3, max.overlaps = 30)+
  geom_text(data= data_prp_complement %>% filter(type=="COMPLEMENT" | secretome_function=="Complement pathway") %>%  count(enrichment),
            aes(x=c(100, 400, 500),y=18, label=paste0("n=",n)))+
  facet_wrap(~enrichment, scales="free_x")+
  theme_classic()
ggsave("outputs/complement.png", width=10, height=3)
```