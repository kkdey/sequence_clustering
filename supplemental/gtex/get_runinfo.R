table_raw = read.delim("SraRunTable.txt")
table = table_raw
table = table[table$Assay_Type_s == "RNA-Seq", ]
attach(table)
table = data.frame(Assay_Type_s, AssemblyName_s, BioSample_s, Experiment_s, 
                   InsertSize_l, LibraryLayout_s, MBases_l, Run_s, SRA_Sample_s, Sample_Name_s, 
                   body_site_s, gap_sample_id_s, gap_subject_id_s)
write.csv(table, "rna_seq_runinfo.csv")
