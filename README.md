# ExpressionPlasticityCanalization
prepare the data
  1. copy the raw count data from Table S3. Save them in a working folder with the codes.
  2. run remove_rRNA_combine_RPM_median_depth_200.pl for flitering expression abundance based on read counts
  3. run total_count_per_cluster.pl to calculate the total splicing events per cluster
  4. run screen_intron_alternative_count.pl to screen for minimal splicing events and combine total events with each event
  5. run combine_cold_warm_condition_count_exp_intron.pl to combine data from different temperature conditions
  
I. plasticity analysis
  1. run simulation_plasticity_null_proportion.R to confirm our approach of identifying consistent plasticity
 i. expresssion abundance
  1. run identify_plastic_gene_expression.pl to identify genes with ancestral plasticity
  2. run identify_PST_outliers_type_expression.pl to identify Pst outliers and classify whether the evolutionary changes are concordant, neutralizing and reversing.
  3. run resample_intron_half_w15_plas_evolv_reidentify_PST.pl to test the patterns of concordant without the potential artifact.
 ii. intron usage
  1. run identify_plastic_intron_usage.pl
  2. run identify_PST_outliers_type_intron.pl to identify Pst outliers and classify whether the evolutionary changes are concordant, neutralizing and reversing.
  3. run resample_intron_half_w15_plas_evolv_reidentify_PST.pl to test the patterns of concordant without the potential artifact.

II. canalization analysis
  1. run expression_variance_inbred_outbred.pl to estimate the expression variance for different types of samples for different genes
  2. run intron_variance_inbred_outbred.pl to estimate the expression variance for different types of samples for different genes
