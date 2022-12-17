
library(tidyverse)
library(dplyr)
library(tidyMB)
library(vegan)
library(broom)
library(DESeq2)
library(biobroom)
library (pander) 




counts <- read_rds("asv_table16ParentsM2novaF.rds")

tax <- read_rds("tax16ParentsMnova2F.rds") %>% 
   as.data.frame() %>% rownames_to_column("Seq") %>% as_tibble() %>% 
   mutate(ASVID = paste0("ASV", 1:nrow(.))) %>% 
   mutate(Phylum2 = ifelse(Phylum == "Proteobacteria", paste("Proteobacteria", as.character(Class)), as.character(Phylum)))
colnames(counts) <- tax$ASVID

total_map <- read_csv("ParentsM_16S_MAP_Nova_F.csv")

# Fixing colnames of counts
#DADA2 uses the ASV sequences as the column names. 
#It takes up a lot of space - so we will convert them to ASV1, ASV2, etc.

tax$ASVID <- paste0("ASV", 1:nrow(tax))
colnames(counts) <- tax$ASVID

## Remove organellar reads from dataframe
organellar_ASV <- tax %>% 
   filter(grepl("Mito", Family) | grepl("Chloropl", Order))
dim(counts)
## Filter out samples with less than 2000 reads
counts2 <- counts[rowSums(counts) > 2000, !colnames(counts)%in%organellar_ASV$ASVID]
dim(counts2)
#Also remove ASVs detected in less than 1% of the samples#
prevalence_function <- function(x){
  return(length(x[x > 0]) / length(x))
}

counts3 <- counts2[,apply(counts2, 2, prevalence_function) > 0.01] 
#dim(counts3)

parentsm_data <- counts3 %>% 
   as.data.frame() %>% 
   rownames_to_column("SampleID") %>% 
   gather(ASVID, value, -SampleID) %>% 
   inner_join(total_map) %>% 
   group_by(SampleID) 
   


## Spread the data into raw count table
## Also remove ASVs detected in less than 5% of the samples

spread_data <- parentsm_data  %>% 
   dplyr::select(SampleID, Compartment, Gen, Location, TRT, ASVID, value) %>% 
   group_by(ASVID) %>% 
   spread(ASVID, value) %>% 
   group_by(SampleID) %>%
   mutate(rank = 1:n()) %>% 
   filter(rank == 1)

## Extract the sample information from the original  dataset
spread_info <- spread_data %>% 
   dplyr::select(SampleID, Compartment, Location, Gen, TRT) %>% 
   mutate(condition = paste(TRT, Compartment, sep = "_")) %>% 
  
   column_to_rownames('SampleID') 
## Make a count table spread into relative abundance

spread_df <- spread_data %>% as_tibble() %>% 
   dplyr::select(-c(Compartment, Location, Gen,  TRT)) %>% 
   column_to_rownames('SampleID')

spread_df_ra <- spread_df / rowSums(spread_df) * 1000

## Field analysis
### Shannon diversity
# This shannon object will be used later for the greenhouse experiments as well
shannon <- as.data.frame(diversity(spread_df)) %>% 
   rownames_to_column("SampleID") %>% 
   dplyr::rename(shannon = 2) %>% 
   inner_join(total_map)

field_diversity <- shannon %>% 
   filter(TRT %in% c("Corpus", "Austin")) %>% 
   mutate(Compartment = fct_relevel(Compartment, "Root", "Rhizo", "Soil"))



plot<-ggplot(field_diversity, aes(Compartment, shannon, color = TRT)) +
   geom_boxplot(size=1.1) +
   geom_point(size =3.5, position = position_jitterdodge(jitter.width = 0.2)) +
   theme_minimal(base_size=40) +
   coord_flip() +
   scale_color_manual(values = c("red", "steelblue"))
plot2<-plot+labs(x="Compartment", y="Shannon Entropy")
plot2


#PCOA




field_info <- spread_info%>% filter(TRT %in% c("Austin", "Corpus")) %>% rownames_to_column("SampleID") 

field_ra <- spread_df_ra %>% rownames_to_column("SampleID") %>% right_join(field_info %>% dplyr::select(SampleID)) %>% column_to_rownames("SampleID")
field_pc <- capscale(log2(field_ra + 1) ~ 1, distance = "bray")
field_pcoa <- field_pc$CA$u[,1:5] %>% as.data.frame() %>% 
   bind_cols(field_info) %>% 
   mutate(Compartment = fct_relevel(Compartment, "Soil", "Rhizo", "Root")) 


summary(field_pc)$cont %>% 
   as.data.frame %>%       
   round(2) %>% .[1:5]%>%  
   pander::pander() 

   ggplot(field_pcoa, aes(MDS1, MDS2, color = TRT, shape = Compartment)) +
   geom_point(size = 8, alpha = 1) +
   geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
   geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
   theme_minimal(base_size=40) +
   scale_color_brewer(palette = "Set1") +
     labs(x = "PCo1 (20%)", y = "PCo2 (16%)")
      
   ### Stats
   Field_info_ns <- spread_info %>% 
      filter(TRT %in% c("Austin", "Corpus") ) %>% 
      rownames_to_column("SampleID") 
      
   Field_ra_ns <- spread_df_ra %>% rownames_to_column("SampleID") %>% right_join(Field_info_ns %>% dplyr::select(SampleID)) %>% column_to_rownames("SampleID")
   
   adonis2(log2(Field_ra_ns + 1) ~  Gen*Compartment, Field_info_ns)
   
   
   ### Phylum Analyses
   ## Make plot for agglomerated counts by phylum
   
field_phy_plot <- parentsm_data %>% 
      filter(SampleID %in% rownames(spread_info))%>% 
      inner_join(tax, by = "ASVID") %>%
      group_by(SampleID) %>% 
     mutate(RA = value / sum(value) * 100) %>% 
      group_by(Phylum2, SampleID, Compartment, TRT) %>% 
      filter(TRT %in% c("Austin", "Corpus"))%>% 
      summarise(phy_tot = sum(RA)) %>% 
      group_by(Phylum2) %>% nest() %>% 
      mutate(phy_all = map_dbl(data, ~sum(.x$phy_tot))) %>% ungroup() %>% 
      top_n(13, phy_all) %>% unnest(data) %>% 
      mutate(Compartment = fct_relevel(Compartment, "Soil", "Rhizo", "Root")) %>% 
      group_by(Compartment, TRT, Phylum2) %>% 
      summarise(mean_tot = mean(phy_tot)) %>% 
      mutate(Phylum2 = fct_relevel(Phylum2, "Verrucomicrobia", after = 8)) %>% 
      mutate(Phylum2 = fct_relevel(Phylum2, "Thaumarchaeota", after = 8)) %>% 
      filter(!is.na(TRT)) 
   
   
   
   
   ggplot(field_phy_plot, aes(Compartment, mean_tot, fill = Phylum2)) + geom_bar(stat = "identity") +
      scale_fill_manual(values = c(RColorBrewer::brewer.pal(12,"Paired")[-c(5,6)], RColorBrewer::brewer.pal(5, "Reds"), "dodgerblue")) +
      facet_grid(.~TRT, scales = "free_x") +
      labs(x = "Sample", y = "Relative abundance (%)") +
      theme_minimal(base_size=40) +
      theme(axis.ticks.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
   
   ## run linear model on phylum relative abundances to ask which phyla are differentially abundant by location within each compartment
  
   field_phy_lm <- parentsm_data %>% 
      filter(SampleID %in% rownames(spread_info)) %>% 
      inner_join(tax, by = "ASVID") %>%
      group_by(SampleID) %>% 
      mutate(RA = value / sum(value) * 1000) %>% 
      filter(Location == "Field") %>% 
      group_by(Phylum2, SampleID, Compartment, Gen, TRT, Location) %>% 
      summarise(phy_tot = sum(RA)) %>% 
      group_by(Phylum2, Compartment) %>% nest() %>% 
      mutate(mod = map(data, ~tidy(lm(log2(phy_tot + 1) ~ TRT, .)))) %>% 
      dplyr::select(Compartment, Phylum2, mod) %>%
      unnest(mod)
   
   ## Make a heatmap of phylum estimates from linear model
   phy_field_lm_plot <- field_phy_lm %>% 
      filter(term != "(Intercept)") %>% 
      ungroup() %>% 
      group_by(Compartment) %>% 
      mutate(p.adj = p.adjust(p.value)) %>% 
      mutate(sig = ifelse(p.adj < 0.05, "sig", NA)) %>% ungroup() %>% 
      mutate(Compartment = fct_relevel(Compartment, "Soil", "Rhizo", "Root")) %>% 
      filter(!is.na(Phylum2)) 
   
      ggplot(phy_field_lm_plot, aes(Compartment, Phylum2, fill = estimate, alpha = sig)) +
      geom_tile(alpha = 1) +
      geom_tile(data = . %>% filter(sig == "sig"), size = 0.5, fill = NA, color = "black") +
      scale_fill_gradient2(low = "red", high = "steelblue") +
      scale_color_manual(values = c("black")) +
         theme_minimal(base_size=30) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      coord_equal()
   field_phy_plot
   phy_field_lm_plot
   

# different aboundance

   dds_soil <- DESeqDataSetFromMatrix(t(spread_df), spread_info, design = ~ condition)
   dds_soil <- DESeq(dds_soil) ## This object will also be used for the greenhouse analysis
   
   
   field_soil <- results(dds_soil, contrast = c("condition", "Austin_Soil", "Corpus_Soil")) %>%tidy() %>% tibble()%>%
      mutate(Compartment = "Soil")
   field_rhizo <- results(dds_soil, contrast = c("condition", "Austin_Rhizo", "Corpus_Rhizo")) %>% tidy() %>% tibble()%>%
      mutate(Compartment = "Rhizo")
   field_root <- results(dds_soil, contrast = c("condition", "Austin_Root", "Corpus_Root")) %>% tidy() %>% tibble()%>%
      mutate(Compartment = "Root")

   field_da_plot <- bind_rows(field_soil, field_rhizo, field_root) %>% 
   filter(p.adjusted < 0.05) %>% 
      inner_join(tax, by = c("gene" = "ASVID")) %>% 
      mutate(direction = ifelse(estimate > 0, "Austin", "Corpus")) %>% 
      dplyr::count(Compartment, Phylum, direction) %>% 
      mutate(n2 = ifelse(direction == "Austin", -n, n)) %>% 
      mutate(Phylum = fct_reorder(Phylum, n),
             Compartment = fct_relevel(Compartment, "Soil", "Rhizo", "Root")) %>% 
      filter(!is.na(Phylum)) %>% 
      group_by(Phylum) %>% filter(sum(n) > 5) 
   
    
      ggplot(field_da_plot, aes(Phylum, n2, fill = direction)) +
      geom_bar(stat = "identity") +
      facet_grid(.~Compartment) +
      scale_fill_manual(values = c("red", "steelblue")) +
      coord_flip() +
      theme_minimal(base_size=40)+
         labs(x = "Phylum", y = "number of enriched ASVs") 
      
      
      
      
      # Greenhouse analysis
      ### PCoA
      gh_info <- spread_info %>% 
         filter(Gen %in% c ("FIL", "HAL") )%>% 
         filter(TRT %in% c ("HI", "MHI", "FI", "MFI") )%>% 
         rownames_to_column("SampleID")
      gh_ra <- spread_df_ra %>% rownames_to_column("SampleID") %>% right_join(gh_info %>% dplyr::select(SampleID)) %>% column_to_rownames("SampleID")
      gh_pc <- capscale(log2(gh_ra + 1) ~ 1, distance = "bray")
      gh_pcoa <- gh_pc$CA$u[,1:5] %>% as.data.frame() %>% 
         bind_cols(gh_info) %>% 
         mutate(Compartment = fct_relevel(Compartment, "Soil", "Rhizo", "Root")) 
      
   
      summary(gh_pc)$cont %>% 
         as.data.frame %>%       
         round(2) %>% .[1:5]%>%
         pander::pander() 
     
      
      
      ggplot(gh_pcoa ,aes(MDS1, MDS2, color = TRT, shape = Compartment)) +
         geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
         geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
         geom_point(size = 8, alpha = 0.8) +
         theme_minimal(base_size=40) +
         scale_color_manual(values = c("#1F78B4","#E31A1C",  "#A6CEE3",  "#FB9A99"  )) +
         labs(x = "PCo1 (21%)", y = "PCo2 (12%)")
      
      
      ### Stats
      gh_info_ns <- spread_info %>% 
         filter(Gen %in% c ("FIL", "HAL") )%>% 
         filter(TRT %in% c("MHI", "HI", "MFI", "FI") ) %>% 
         rownames_to_column("SampleID") 
      
      gh_ra_ns <- spread_df_ra %>% rownames_to_column("SampleID") %>% right_join(gh_info_ns %>% dplyr::select(SampleID)) %>% column_to_rownames("SampleID")
      
      adonis2(log2(gh_ra_ns + 1) ~  Gen*TRT*Compartment, gh_info_ns)
      
      ### Test for genotype effect without soil
      gh_info_ns <- spread_info %>% 
         filter(Gen %in% c ("FIL", "HAL") )%>% 
         filter(TRT %in% c("MHI", "HI", "MFI", "FI") ) %>% 
         rownames_to_column("SampleID") %>% 
         filter(Compartment != "Soil")
      gh_ra_ns <- spread_df_ra %>% rownames_to_column("SampleID") %>% right_join(gh_info_ns %>% dplyr::select(SampleID)) %>% column_to_rownames("SampleID")
      
      adonis2(log2(gh_ra_ns + 1) ~ TRT + Gen + Compartment, gh_info_ns)
      
      
      ### Compare BC dissimilarities between treatment levels
      
      
   df_dist <- vegdist(log2(spread_df + 1)) %>% as.matrix()
   df_dist[upper.tri(df_dist, diag = T)] <- NA
      df_dist_long <- df_dist %>%
         as.data.frame() %>% 
         rownames_to_column("Sample1") %>% 
         gather(Sample2, dist, -Sample1) %>%
         na.omit() %>% 
         inner_join(spread_info %>% rownames_to_column("Sample1"), by = "Sample1") %>% 
         inner_join(spread_info %>% rownames_to_column("Sample2"), by = "Sample2")
      
      gh_dist_comp <- df_dist_long %>% 
         filter(Compartment.x == Compartment.y) %>% 
         mutate(Compartment = Compartment.x) %>% dplyr::select(-c(Compartment.x, Compartment.y)) %>%
         filter(TRT.x%in%c("MHI", "HI", "MFI", "FI") & TRT.y%in%c("MHI", "HI", "MFI", "FI")) %>% 
         filter(TRT.x != TRT.y) %>% 
         mutate(type.x = ifelse(TRT.x%in%c("MHI", "MFI"), "MHI vs MFI", "HI vs FI"),
                type.y = ifelse(TRT.y%in%c("MHI", "MFI"), "MHI vs MFI", "HI vs FI")) %>% 
         filter(type.x == type.y) %>% 
         mutate(Compartment = fct_relevel(Compartment, "Soil", "Rhizo", "Root"))
      
      
         ggplot( gh_dist_comp,aes(Compartment, dist, color = type.x)) +
         geom_boxplot(outlier.alpha = 0) +
         geom_point(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.5, size = 4) +
         scale_color_manual(values = c("black", "gray" )) +
         labs(x = "", y = "Bray-Curtis Dissimilarity") +
         theme_minimal(base_size=40)
      
         
         #Shannon diversity
      
         gh_diversity  <- shannon %>% 
            filter(TRT %in% c("MFI", "FI", "MHI", "HI")) %>% 
            mutate(Compartment = fct_relevel(Compartment, "Soil", "Rhizo", "Root"))
         
         
         
        plot<- ggplot(gh_diversity, aes(Compartment, shannon, color = TRT)) +
            geom_boxplot(size=1.1, outlier.alpha = 0) +
            geom_point(size =2, position = position_jitterdodge(jitter.width = 0.1), size = 1) +
            theme_minimal(base_size=40) +
            coord_flip()  +
            scale_color_manual(values = c("#1F78B4", "#E31A1C", "#A6CEE3", "#FB9A99" ))    
         plot2<-plot+labs(x="Compartment", y="Shannon Entropy")
         plot2
         
         ### Phylum level statistics
      
         gh_phy_plot <- parentsm_data %>% 
            filter(SampleID %in% rownames(spread_info)) %>% 
            inner_join(tax, by = "ASVID") %>%
            group_by(SampleID) %>% 
            mutate(RA = value / sum(value) * 100) %>% 
            group_by(Phylum2, SampleID, Compartment, TRT) %>% 
            summarise(phy_tot = sum(RA)) %>% 
            group_by(Phylum2) %>% nest() %>% 
            mutate(phy_all = map_dbl(data, ~sum(.x$phy_tot))) %>% ungroup() %>% 
            top_n(13, phy_all) %>% unnest(data) %>% 
            mutate(Compartment = fct_relevel(Compartment, "Soil", "Rhizo", "Root")) %>% 
            group_by(Compartment, TRT, Phylum2) %>% 
            summarise(mean_tot = mean(phy_tot)) %>% 
            mutate(Phylum2 = fct_relevel(Phylum2, "Verrucomicrobia", after = 9)) %>% 
            filter(!is.na(TRT) & !TRT%in%c("Corpus", "Austin")) 
         
         
            ggplot(gh_phy_plot,aes(Compartment, mean_tot, fill = Phylum2)) + geom_bar(stat = "identity") +
            scale_fill_manual(values = c(RColorBrewer::brewer.pal(12,"Paired")[-c(5,6)], RColorBrewer::brewer.pal(5, "Reds"), "dodgerblue")) +
            facet_grid(.~TRT, scales = "free_x") +
            labs(x = "Sample", y = "Relative abundance (%)") +
            theme_minimal(base_size=40) +
            theme(axis.ticks.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1))
         
            gh_phy_lm <- parentsm_data %>% 
               filter(SampleID %in% rownames(spread_info)) %>% 
               inner_join(tax, by = "ASVID") %>%
               group_by(SampleID) %>% 
               mutate(RA = value / sum(value) * 1000) %>% 
               group_by(Phylum2, SampleID, Compartment, TRT, Location, Gen) %>% 
               summarise(phy_tot = sum(RA)) %>% 
               filter(!TRT %in% c("Austin", "Corpus")) %>% 
               mutate(type = ifelse(TRT%in%c("HI", "FI"), "treated", "null"),
                      source = ifelse(TRT%in%c("MHI", "HI"), "Austin", "Corpus")) %>%  
               group_by(Phylum2, Compartment, type) %>% nest() %>% 
               mutate(mod = map(data, ~tidy(lm(log2(phy_tot + 1) ~ source, .)))) %>% 
               dplyr::select(Compartment, Phylum2, type, mod) %>%
               unnest(mod)
         
            gh_phy_model_plot <- gh_phy_lm %>% 
               filter(term != "(Intercept)") %>% 
               ungroup() %>% 
               group_by(Compartment) %>% 
               mutate(p.adj = p.adjust(p.value)) %>% 
               mutate(sig = ifelse(p.adj < 0.05, "sig", NA)) %>% ungroup() %>% 
               mutate(Compartment = fct_relevel(Compartment, "Root", "Rhizo", "Soil")) %>% 
               filter(!is.na(Phylum2)) %>% 
               mutate(Phylum2 = gsub("Proteobacteria ", "", Phylum2)) 
            
             
               ggplot( gh_phy_model_plot, aes(Phylum2, Compartment, fill = estimate, alpha = sig)) +
               geom_tile(alpha = 1) +
               geom_tile(data = . %>% filter(sig == "sig"), size = 0.5, fill = NA, color = "black") +
               scale_fill_gradient2(low = "red", high = "steelblue") +
               scale_color_manual(values = c("black")) +
               facet_grid(type~.) +
               theme_minimal(base_size=20) +
               theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
               coord_equal()
            
            gh_phy_plot
            gh_phy_model_plot
         
         # different aboundance
         
         trt_soil <- results(dds_soil, contrast = c("condition", "HI_Soil", "FI_Soil")) %>% tidy() %>% 
            mutate(Compartment = "Soil", type = "TRT")
         trt_rhizo <- results(dds_soil, contrast = c("condition", "HI_Rhizo", "FI_Rhizo")) %>% tidy() %>% 
            mutate(Compartment = "Rhizo", type = "TRT")
         trt_root <- results(dds_soil, contrast = c("condition", "HI_Root", "FI_Root")) %>% tidy() %>% 
            mutate(Compartment = "Root", type = "TRT")
         
         mock_soil <- results(dds_soil, contrast = c("condition", "MHI_Soil", "MFI_Soil")) %>% tidy() %>% 
            mutate(Compartment = "Soil", type = "Mock")
         mock_rhizo <- results(dds_soil, contrast = c("condition", "MHI_Rhizo", "MFI_Rhizo")) %>% tidy() %>% 
            mutate(Compartment = "Rhizo", type = "Mock")
         mock_root <- results(dds_soil, contrast = c("condition", "MHI_Root", "MFI_Root")) %>% tidy() %>% 
            mutate(Compartment = "Root", type = "Mock")
         
         gh_da_asv <- bind_rows(
            trt_soil,
            trt_rhizo,
            trt_root,
            mock_soil,
            mock_rhizo,
            mock_root
         ) %>% 
            dplyr::rename(ASVID = gene) %>% 
            filter(p.adjusted < 0.05) %>% 
            inner_join(tax) %>% 
            mutate(direction = ifelse(estimate > 0, "HI", "MFI")) %>% 
            dplyr::count(Phylum, Compartment, type, direction) %>% 
            mutate(n2 = ifelse(type == "Mock", -n, n)) %>% 
            mutate(Phylum = fct_reorder(Phylum, n),
                   comparison = paste(direction, type)) %>% 
            group_by(Phylum, Compartment) %>% 
            filter(sum(n) > 10) %>% ungroup() %>% 
            mutate(Compartment = fct_relevel(Compartment, "Soil", "Rhizo", "Root")) 
         
         
            ggplot(gh_da_asv, aes(Phylum, n2, fill = comparison)) +
            geom_bar(stat = "identity") + 
            geom_hline(yintercept = 0, size = 0.1) +
            facet_grid(.~Compartment) +
            theme_minimal(base_size=22) +
            labs(x = "Phylum", y = "number of enriched ASVs") +
            coord_flip() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            scale_fill_manual(values = c("#FB9A99", "#E31A1C", "#A6CEE3", "#1F78B4"))
         
         
            
   # Combined Analysis
    ### PCoA
            
           
            combo_pc <- capscale(log2(spread_df_ra + 1) ~ 1)
            head(combo_pc$CA$eig / sum(combo_pc$CA$eig))
            combo_pc_plot <- combo_pc$CA$u[,1:5] %>% as.data.frame() %>% 
               bind_cols(spread_info)#%>% 
           
               summary(combo_pc)$cont %>% 
               as.data.frame %>%    .[1:5]    %>%  
               round(2) %>%
               pander::pander() 
            
            
            ggplot(combo_pc_plot ,aes(MDS1, MDS2, color = TRT, shape = Location)) +
               geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
               geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
               geom_point(size = 8, alpha = 0.8) +
               scale_shape_manual(values = c(1, 16)) +
               theme_minimal(base_size=40) +
               scale_color_manual(values = c("#E31A1C","#1F78B4" ,"#1F78B4" , "#E31A1C",  "#A6CEE3", "#FB9A99")) +
               labs(x = "PCo1 (26%)", y = "PCo2 (9%)")
            
            
            ### Phylum Stats
            
            field_gh_lm <- parentsm_data %>% 
               filter(SampleID %in% rownames(spread_info)) %>% 
               inner_join(tax, by = "ASVID") %>%
               group_by(SampleID) %>% 
               mutate(RA = value / sum(value) * 1000) %>% 
               filter(!TRT%in%c("MFI", "MHI")) %>% 
               group_by(Phylum2, SampleID, Compartment, TRT, Location) %>% 
               summarise(phy_tot = sum(RA)) %>% 
                             group_by(Phylum2, Compartment) %>% nest() %>% 
               mutate(mod = map(data, ~tidy(lm(log2(phy_tot + 1) ~ Location, .)))) %>% 
               dplyr::select(Compartment, Phylum2, mod) %>%
               unnest(mod)
            
            field_gh_comp <- field_gh_lm %>% 
               filter(term != "(Intercept)") %>% 
               ungroup() %>% 
               group_by(Compartment) %>% 
               mutate(p.adj = p.adjust(p.value)) %>% 
               mutate(sig = ifelse(p.adj < 0.05, "sig", NA)) %>% ungroup() %>% 
               mutate(Compartment = fct_relevel(Compartment, "Soil", "Rhizo", "Root")) %>% 
               filter(!is.na(Phylum2)) 
              
               
                ggplot( field_gh_comp, aes(Phylum2, Compartment, fill = estimate, alpha = sig)) +
               geom_tile(alpha = 1) +
                   geom_tile(data = . %>% filter(sig == "sig"), size = 0.5, fill = NA, color = "black") +
               scale_fill_gradient2(low = "gold", high = "darkorchid2") +
               scale_color_manual(values = c("black")) +
               scale_alpha_manual(values = c(0,1)) +
              
                   theme_minimal(base_size=25) +
               theme(axis.text.x = element_text(angle = 45, hjust = 1))+
               labs(x = "Phylum", y = "") +
               coord_equal()


            
    ## similarity in field and greenhouse
                
                parentsm_data2 <- counts3 %>% 
                  as.data.frame() %>% 
                  rownames_to_column("SampleID") %>% 
                  gather(ASVID, value, -SampleID) %>% 
                  inner_join(total_map) %>% 
                  group_by(SampleID) %>% 
                mutate(RA = value / sum(value) * 100)
                
       # Rhizosphere         
               
                f_gh_mod_Rhizo <- parentsm_data2 %>%
                  filter(Compartment%in%c("Rhizo")) %>%
                  mutate(TRT2 = ifelse(TRT%in%c("Austin", "HI"), "Austin", "Corpus")) %>%
                  filter(!TRT%in%c("MFI", "MHI")) #%>%
                  group_by( Location, ASVID) %>%
                 filter(sum(RA > 0)/ n() > 0.3) %>%
                  nest() %>%
                  mutate(mod = map(data, ~tidy(lm(log2(RA + 1) ~ TRT2, .)))) %>%
                  unnest(mod) %>% select(-data)
                
    ## Count differentially abundant microbes between AUS/Corpus in field and greenhouse inoculation
                
                f_gh_mod_Rhizo %>%
                  ungroup() %>%
                  filter(term != "(Intercept)") %>%
                  group_by(Location) %>%
                  mutate(padj = p.adjust(p.value, "BH")) %>%
                  mutate(direction = ifelse(estimate < 0, "Austin", "Corpus")) %>%
                  filter(padj < 0.05) %>%
                  ungroup() %>%
                  dplyr::count(Location, direction)
          
                
    ## Find overlapping in directionality between the greenhouse and field settings
                f_gh_mod_RhizoF <-   f_gh_mod_Rhizo %>%
                  ungroup() %>%
                  filter(term != "(Intercept)") %>%
                  group_by(Location) %>%
                  mutate(padj = p.adjust(p.value, "BH")) %>%
                  mutate(direction = ifelse(estimate < 0, "Austin", "Corpus")) %>%
                  filter(padj < 0.05) %>%
                  ungroup() %>%
                  dplyr::count(ASVID, direction) %>%
                  filter(n == 2)  %>%            
                 print (n=39)
         
             
             
    # Root
        
             
             f_gh_mod_Root <- parentsm_data2 %>%
               filter(Compartment%in%c("Root")) %>%
              
               mutate(TRT2 = ifelse(TRT%in%c("Austin", "HI"), "Austin", "Corpus")) %>%
               filter(!TRT%in%c("MFI", "MHI")) %>%
              
               group_by(Location, ASVID) %>%
               filter(sum(RA > 0)/ n() > 0.3) %>%
               nest() %>%
               mutate(mod = map(data, ~tidy(lm(log2(RA + 1) ~ TRT2, .)))) %>%
               unnest(mod) %>% select(-data)
             
             ## Count differentially abundant microbes between AUS/Corpus in field and greenhouse inoculation
             
             f_gh_mod_Root %>%
               ungroup() %>%
               filter(term != "(Intercept)") %>%
               group_by(Location) %>%
               mutate(padj = p.adjust(p.value, "BH")) %>%
               mutate(direction = ifelse(estimate < 0, "Austin", "Corpus")) %>%
               filter(padj < 0.05) %>%
               ungroup() %>%
               dplyr::count(Location, direction)
             
             
             ## Find overlapping in directionality between the greenhouse and field settings
             f_gh_mod_Rootf <-        f_gh_mod_Root %>%
               ungroup() %>%
               filter(term != "(Intercept)") %>%
               group_by(Location) %>%
               mutate(padj = p.adjust(p.value, "BH")) %>%
               mutate(direction = ifelse(estimate < 0, "Austin", "Corpus")) %>%
               filter(padj < 0.05) %>%
               ungroup() %>%
               dplyr::count(ASVID, direction) %>%
               filter(n == 2)
            
             
  
             
     # Soil
             
             
             f_gh_mod_Soil <- parentsm_data2 %>%
               filter(Compartment%in%c("Soil")) %>%
               
               mutate(TRT2 = ifelse(TRT%in%c("Austin", "HI"), "Austin", "Corpus")) %>%
               filter(!TRT%in%c("MFI", "MHI")) %>%
             
               group_by(Location, ASVID) %>%
              filter(sum(RA > 0)/ n() > 0.3) %>%
               nest() %>%
               mutate(mod = map(data, ~tidy(lm(log2(RA + 1) ~ TRT2, .)))) %>%
               unnest(mod) %>% select(-data)
             
             ## Count differentially abundant microbes between AUS/Corpus in field and greenhouse inoculation
             
             f_gh_mod_Soil %>%
               ungroup() %>%
               filter(term != "(Intercept)") %>%
               group_by(Location) %>%
               mutate(padj = p.adjust(p.value, "BH")) %>%
               mutate(direction = ifelse(estimate < 0, "Austin", "Corpus")) %>%
               filter(padj < 0.05) %>%
               ungroup() %>%
               dplyr::count(Location, direction)
             
             
             ## Find overlapping in directionality between the greenhouse and field settings
             f_gh_mod_Soil<-   f_gh_mod_Soil %>%
               ungroup() %>%
               filter(term != "(Intercept)") %>%
               group_by(Location) %>%
               mutate(padj = p.adjust(p.value, "BH")) %>%
               mutate(direction = ifelse(estimate < 0, "Austin", "Corpus")) %>%
               filter(padj < 0.05) %>%
               ungroup() %>%
               dplyr::count(ASVID, direction) %>%
               filter(n == 2)  %>%            
               print (n=39)        
            
                
               
            
            
            
            
         

