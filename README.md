#bcr_in_three_cancers
This repository contains data and code used for analysis and visualisations in Krasik et al (2024) Systematic evaluation of intratumoral and peripheral BCR repertoires in three cancers eLife 13:RP89506 https://doi.org/10.7554/eLife.89506.2

Please be aware that the code is given as a working example for a specific set of parameters and paths. Different parameters might have been used to produce similar plots, e.g. for another patient or for another tissue etc.

The following are the code used to produce figures and the corresponding data:

Fig 1B replicates Scatterplots.py, BCR_CDR3_data/5RACE/All_data_filtered/

Fig 1B fragments Scatterplots.py, BCR_CDR3_data/5RACE/Pooled_Replicas/mixcr/

Fig 2A Pipelines/Correlation_plot_attempt.py, BCR_CDR3_data/5RACE/Pooled_Samples_aa/mixcr/

Fig 2B bulk_vs_PCsort.py, PwDistances/ (calculated from BCR_CDR3_data/5RACE/Top_for_patient/Pooled_Replicas/VDJ/)

Fig 2C bulk_vs_PCsort.py, BCR_CDR3_data/5RACE/IsotypeComposition/

Fig 3A Overlap_comparison_all_isotypes.py, PwDistances (calculated from BCR_CDR3_data/5RACE/Pooled_Samples_aa/VDJ/)

Fig 3C Norm_SW_clonality.py, DiversityStats (calculated from BCR_CDR3_data/5RACE/All_data_filtered/)

Fig 3D Norm_SW_clonality.py, DiversityStats (calculated from BCR_CDR3_data/5RACE/All_data_filtered/)

Fig 3E Overlap_by_isotypes.py, PwDistances (calculated from  /home/sofyakrasik/BCR_CDR3_data/5RACE/Pooled_Samples_aa/VDJ/)

Fig 3F Overlap_by_isotypes.py, PwDistances (calculated from  /home/sofyakrasik/BCR_CDR3_data/5RACE/Pooled_Samples_aa/VDJ/)

Fig 3G Isotype_correlation.py, BCR_CDR3_data/5RACE/All_data_filtered/

Fig 3H Isotype_correlation.py, BCR_CDR3_data/5RACE/All_data_filtered/

Fig 4A CDR3_length_comparison.py, BCR_CDR3_data/5RACE/Pooled_Replicas_aa/VDJ/

Fig 4B CDR3_length_comparison.py, BCR_CDR3_data/5RACE/Pooled_Replicas_aa/VDJ/

Fig 4C AA_R_heatmap.py, BCR_CDR3_data/5RACE/AminoAcidPropertiesStats/Pooled_Replicas/

Fig 4D AA_R_heatmap.py, BCR_CDR3_data/5RACE/AminoAcidPropertiesStats/Pooled_Replicas/

Fig 4E Number_of_SHM.py, Clonal_lineages/ChangeOut/

Fig 5B Cancer_Triangle_Plot.py, Clonal_lineages/ChangeOut/

Fig 5C Cancer_Triangle_Plot.py, Clonal_lineages/ChangeOut/

Fig 5D Cancer_Triangle_Plot.py, Clonal_lineages/ChangeOut/

Fig 6A Heterogeneity_quantified.py, PwDistances/SeparateReplicas/

Fig 6C Heterogeneity_triangles_new.py, Clonal_lineages/ChangeOut/

Fig 6D Pipeline_Working.py (plotting features may vary a lot on size, shape, colour etc.), metadata example is in Clonal_lineages/Metadata/

Fig 7A Heterogeneity_quantified.py, PwDistances/SeparateReplicas_with_singles/

Fig 7B Heterogeneity_triangles_new.py, Clonal_lineages/ChangeOut/ (CDR3 with singles)

Fig 7D Heterogeneity_triangles_new_with_isotypes.py, Clonal_lineages/ChangeOut/ (CDR3 with singles)

Fig 8A Expanded_clonotypes_Segment_level.py - to identify expanded clonotypes with EdgeR, result in BCR_CDR3_data/5RACE/Expanded_clonotypes/segment_level (calculated from BCR_CDR3_data/5RACE/Top_for_patient/Separate_Replicas)

Fig 8B Expanded_clonotypes_Segment_level.py - to identify expanded clonotypes with EdgeR, result in BCR_CDR3_data/5RACE/Expanded_clonotypes/segment_level (calculated from BCR_CDR3_data/5RACE/Top_for_patient/Separate_Replicas)

Fig 8D Pipeline_Working.py (plotting features may vary a lot on size, shape, colour etc.), metadata example is in Clonal_lineages/Metadata/

Fig 8E Pipeline_Working.py (plotting features may vary a lot on size, shape, colour etc.), metadata example is in Clonal_lineages/Metadata/

Fig 8F Expanded_clonotypes_Tissue_Level_with_Alternative_pooling.py  - to identify expanded clonotypes with EdgeR, results in BCR_CDR3_data/5RACE/Expanded_clonotypes_Alternative (calculated from BCR_CDR3_data/5RACE/Alternative_Replicas_Pooling_aa/)



Clonotype sharing visualisations were created manually using Cytoscape.

For manipulations on original data (e.g. isotype proportions, sample integration, selection of top clonotypes etc.), see Core_folders_creator.py

