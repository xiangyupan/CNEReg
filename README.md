# CNEReg    
CNEReg is an evolutionary Conserved Non-coding Element interpretation method by gene Regulatory network to integrate paired time-course expression and chromatin accessibility data with the available public data, including comparative genomics data and comparative transcriptomic data. CNEReg aims to systematically fill the gap between conserved non-coding elements (CNEs) and its significantly impacted morphology in evolution. This is done by reconstructing a developmental regulatory network by paired time series of paired gene expression and chromatin accessibility data.Particularly in sheep CNEs are RSCNEs and morphology is the innovation of rumen organ, which is further denoted by the set of rumen specific genes RSEGs. We reconstruct gene regulatory network during rumen development to systematically understand how the TFs regulate genes via batteries of RSCNEs, which over development, lead to the cell type specific activation of RSEGs. The main idea of CNEReg is to focus on those toolkit TFs as major players in evo-devo studies to study how those TFs are regulated by RSCNEs and how they utilize RSCNEs to regulate RSEGs. CNEReg models the expression of target genes (TG) conditional on chromatin accessibility of RSCNEs and expression of transcription factors (TF).  

____________________________________      
# 1. Data availablility    
# (1) Raw data    
The raw data files are the original files containing reads and quality scores, as generated by the Illumina sequencing instrument. Raw Data File Format if FASTQ. All the RAW data and processed data produced in this study have been permanently deposited at the Sequence Read Archive (SRA) under project number PRJNA485657. The useres and reviewers can access the raw data with fastq files via URL: <https://submit.ncbi.nlm.nih.gov/subs/sra/SUB7461604/>.  
# (2) Processed data      
CNEReg model requires input as sample matched time-series RNA-seq and ATAC-seq, and RSCNEs with conservation scores from public data.       
The final processed data are defined as the data on which the conclusions in the related manuscript are based.    
  1. Peak files with quantitative openness data with a format bed and txt files for ATAC-Seq data.    
  2. The normalized gene expression profile output from Stringtie for RNA-seq data.   
  3. The bed file of RSCNEs which consist of chromosome, start and end coordinates with the corresponding conservation score.
# 2. Processing data    
   

# (1) The "heatmap&PCA_of_Chromatin_accessibility.R" script was used to draw the correlation heatmap and PCA for chromatin accessibility.   
The input file "openness_RUMEN.csv" includes the normalized openness matrix of each active-RSCNEs of each time point.    
`Rscript heatmap&PCA_of_Chromatin_accessibility.R openness_RUMEN.csv heatmap_of_Chromatin_accessibility.pdf PCA_of_Chromatin_accessibility.pdf`   
# (2) The "heatmap&PCA of Gene Expression.R" script was used to draw the correlation heatmap and PCA for gene expression.   
The first input file "RNA-seq_RUMEN.csv" includes the gene FPKM matrix file of rumen and esophagus during development.    
The second input file "sif.xlsx" stored the batch effect info.   
Two output files are in PDF format.   
`Rscript heatmap&PCA of Gene Expression.R RNA-seq_RUMEN.csv sif.xlsx heatmap_of_Gene_expression.pdf  PCA_of_Gene_expression.pdf`   
# (3) The "JMscore.R" script was used to calculate the JSD scores and median expression data of all the genes in 50 tissues and combined them into JMscore to measure the specificity of genes in different tissues with the subset of JMscore of 18 TTF.   
The input file "sheep832.csv" includes the gene expression level of each gene in each sample from 830 sheep samples.   
The first output file "JSD_sheep.csv" includes the Jensen-Shannon divergence of ench gene in each tissue.   
The second output file "median_sheep.csv" includes the median expression level of ench gene in each tissue.   
The third output file "JMscore_sheep.csv" includes the JMscore of ench gene in each tissue.   
`Rscript JMscore.R sheep832.csv JSD_sheep.csv JSD_sheep.csv JMscore_sheep.csv`    
# (4) The "18TTF-JMscore-cluster.R" script was used to darw the hierarchical clustering diagram.      
The first input file "18TTFs-JSD.csv" includes the Jensen-Shannon divergence score of 18 TTFs in each tissue.   
The second input file "18TTFs-median.csv" includes the median expression level of 18 TTFs in each tissue.   
The output file "Phylogeny_of_50_tissues.pdf" is the phylogeny tree of 50 tissues clustered by 18 TTFs.     
# (5) The "upstream(lasso).R" script was used to generate the TFs and RSCNEs under the linear regression model.        
The first input file "RNA-seq-TF.csv" includes the expression level of each TF at each time point.      
The second input file "openness.csv" includes the openness of each RSCNE at each time point.      
The third input file "19TTF-upbinding.csv" includes the TF binding information of each TTF which predicted by homer.      
The fourth input file "19TTF-RSCNE(LASSO).csv" is the RSCNEs selected by LASSO.   
The outputfile "19TTF-up-lm.csv" includes the TFs and RSCNEs for linear regression.   
# (6) The "upstream(lm).R" scrpit was used to generate the RSCNEs of 19TTFs for the upstream network.     
The first input file "RNA-seq-TF.csv" includes the expression level of each TF at each time point.        
The second input file "openness.csv" includes the openness of each RSCNE at each time point.        
The third input file "19TTF-up-lm(cor0.6).csv" includes the TFs selected by cor(TF,TTF)>0.6.    
The otput file "19TTF-up-RSCNE.csv" includes the RSCNEs remained in the upstream network.   
# (7) The "downstream.R" script was used to calculate the TTFs downstream network.    
The first input file "openness.csv" includes the openness of each RSCNE at each time point.   
The second input file "RNA-seq.csv" includes the expression level of each gene at each time point.    
The third input file "17TTF-downstream.csv" includes the regulatory relationship of TTFs predicted by homer.    
The output file "downstream_network.csv" includes the regulatory network of TTFs in the downstream.   
# (8) The "functional influence_upstream.R" script was used to calculate the functional influence of active-RSCNEs in the TTF upstream network.     
The input file "ATAC_average.xlsx" is the accessibility level of each active-RSCNE at each time point.
The input file "RNA_average.xlsx"  is the gene expression level of each TF at each time point.    
The input file "Binding&Correlation_up.csv" is the motif bingding strength ![image](https://github.com/xiangyupan/CNEReg/blob/main/CodeCogsEqn3.svg) and spearman correlation of ![image](https://github.com/xiangyupan/CNEReg/blob/main/CodeCogsEqn1.svg) and ![image](https://github.com/xiangyupan/CNEReg/blob/main/CodeCogsEqn2.svg)    
The input file "URSCNE.xlsx" is the active-RSCNEs in the upstream network.    
The input file "type1.conservativeS.txt" is the conservation scores of type I active-RSCNES.    
The input file "type2.conservativeS.txt" is the conservation scores of type II active-RSCNES.   
The output file "up-FI.csv" is the functional influence of active-RSNCEs in the TTF upstream network.   
`Rscript functional_influence_upstream.R ATAC_average.xlsx RNA_average.xlsx Binding&Correlation_up.csv URSCNE.xlsx type1.conservativeS.txt type2-a.conservativeS.txt up-FI.csv`     
# (9) The "functional influence_downstream.R" script was used to calculate the functional influence of active-RSCNEs in the TTF downstream network.   
The input file "ATAC_average.xlsx" is the accessibility level of each active-RSCNE at each time point.    
The input file "RNA_average.xlsx"  is the gene expression level of each TF at each time point.        
The input file "Binding&Correlation_up.csv" is the motif bingding strength ![image](https://github.com/xiangyupan/CNEReg/blob/main/CodeCogsEqn3.svg) and spearman correlation of ![image](https://github.com/xiangyupan/CNEReg/blob/main/CodeCogsEqn1.svg) and ![image](https://github.com/xiangyupan/CNEReg/blob/main/CodeCogsEqn2.svg)    
The input file "DRSCNE.xlsx" is the active-RSCNEs in the downstream network.    
The input file "type1.conservativeS.txt" is the conservation scores of type I active-RSCNES.    
The input file "type2.conservativeS.txt" is the conservation scores of type II active-RSCNES.   
The output file "down-FI.csv" is the functional influence of active-RSNCEs in the TTF downstream network.     
`Rscript functional_influence_downstream.R ATAC_average.xlsx RNA_average.xlsx Binding&Correlation_up.csv DRSCNE.xlsx type1.conservativeS.txt type2-a.conservativeS.txt down-FI.csv`     
# (10) The "functional influence_diffNetwork.R" script was used to calculate the functional influence of active-RSCNEs in the differential subnetwork between rumen and esophagus.   
The input file "diff-single-system.txt" is the regulatory strength of each TF-RSCNE-TG pair at each time point.     
The input file "DiffRSCNE.xlsx" stored active-RSCNEs in the differential subnetwork between rumen and esophagus.      
The input file "type1.conservativeS.txt" is the conservation scores of type-I active-RSCNEs.      
The input file "type2-a.conservativeS.txt" is the conservation scores of type-II active-RSCNEs.   
The output file "diff-FI.csv" is the functional influence of active-RSCNEs in the differential subnetwork between rumen and esophagus.      
`Rscript functional influence_diffNetwork.R diff-single-system.txt DiffRSCNE.xlsx type1.conservativeS.txt type2-a.conservativeS.txt diff-FI.csv`    
  
#  Reference    
Xiangyu Pan, Zhaoxia Ma, Xinqi Sun, Hui Li, Tingting Zhang, Zhao Chen, Nini Wang, Wing Hung Wong, Wen Wang, Yu Jiang, Yong Wang. Interpreting ruminant specific conserved non-coding elements by reconstructing developmental regulatory network. (In submission).   

* We'd love to hear from you. If you have any questions, please don't be hestitate to contact the author of this manuscript: [pan_xiangyu@nwafu.edu.cn](pan_xiangyu@nwafu.edu.cn)   
