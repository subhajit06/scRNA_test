#### sample information
library(data.table)
library(dplyr)
library(stringr)

df_id_0 = fread("../ROSMAP_scRNA_data/Metadata/id_mapping.csv")
L = dim(df_id_0)[1]

df_id = NULL
df_id$projid = unique(df_id_0$projid[seq(1, L, 12)])
df_id$Subject = unique(df_id_0$Subject[seq(1, L, 12)])
df_id = as.data.frame(df_id)
df_id$Subject = as.character(df_id$Subject)

ros_clinical_df_1 = fread("../../../ROSMAP/nature_paper_48_subj_AD_info.csv")
ros_clinical_df_2 = fread("../../../ROSMAP/nature_paper_48_subj_AD_info_pathology.csv")
ros_clinical_df_3 = fread("../../../ROSMAP/nature_paper_48_subj_ind_id.csv")

ros_clinical_df_1_2 = inner_join(ros_clinical_df_1, ros_clinical_df_2, by = "Subject") 
saveRDS(ros_clinical_df_1_2, file="./ros_clinical_df.Rds")

ros_clinical_df = ros_clinical_df_1_2 %>%
                      select(Subject, msex.y, pathologic_diagnosis_of_AD, amyloid.y)

df_ros_clinical_projid = inner_join(ros_clinical_df, df_id, by = "Subject")

tmp1 = fread("../ROSMAP_scRNA_data/Genotypes/Raw/ROSMAP_arrayGenotype.csv")
tmp2 = NULL
tmp2$projid = unlist(lapply(tmp1$gwas_id, function(x){substr(x, 4, str_length(x))})) 
tmp2$genotype_present = rep(1, length(tmp2$projid))
tmp2 = as.data.frame(tmp2)
tmp2$projid = as.integer(as.character(tmp2$projid))
df_ros_with_genotype = left_join(df_ros_clinical_projid, tmp2, by = "projid")
df_ros_with_genotype$genotype_present[is.na(df_ros_with_genotype$genotype_present) == TRUE] = 0

table(df_ros_with_genotype$pathologic_diagnosis_of_AD, df_ros_with_genotype$genotype_present)
### 18 out of those 36 ppl have AD+ 
#0  1
#NO   6 18   
#YES  6 18

saveRDS(df_ros_with_genotype, "./df_ros_with_genotype.Rds")

keep_ind = NULL
geno_fam_df = fread("../ROSMAP_scRNA_data/Genotypes/Raw/ROSMAP_arrayGenotype.fam",header = FALSE)
within_fam_id = as.list((df_ros_with_genotype %>% mutate(ROS_id = sprintf("ROS%d", projid)))$ROS_id)

sel_gegno_sample_having_scRNA = geno_fam_df %>% filter(V2 %in% within_fam_id)

keep_ind$fam_id = sel_gegno_sample_having_scRNA$V1
keep_ind$within_fam_id = sel_gegno_sample_having_scRNA$V2
keep_ind = as.data.frame(keep_ind)
fwrite(keep_ind, "./keep_ind_id.txt",quote=FALSE,sep="\t", col.names = FALSE)
### plink --bfile ROSMAP_arrayGenotype --keep ../../../rosmap_code/keep_ind_id.txt --make-bed --out data_keep_scRNA
### 750173 variants and 36 people pass filters and QC


