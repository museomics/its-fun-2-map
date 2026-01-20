its_outcome <- function(df){

	df$decision_description <- NA
	df$Final_outcome <- NA
	df$Final_contig_desc <- NA
	df$Final_contig <- NA
	
	### Scenario 1,3 & 4
	s1 <- which(df$blast_its2_correct_taxonomy	== "PASS" & 
				df$blast_its1_correct_taxonomy == "PASS" &
				df$overlapping_coords == "NO" &
				df$same_contig == "YES" &
				df$extraction_ITS_complete == "NO"
	)
	
	
	df[s1, "decision_description"] <- "Putative ITS region"
	df[s1, "Final_contig"] <- df$blast_its2_contig_path[s1]
	df[s1, "Final_contig_desc"] <- "putative_ITS_region"
	df[s1, "Final_outcome"] <- "PASS"
	
	
	### Scenario 2
	s2 <- which(df$blast_its2_correct_taxonomy	== "PASS" & 
				df$blast_its1_correct_taxonomy == "PASS" &
				df$overlapping_coords == "NO" &
				df$same_contig == "YES" &
				df$extraction_ITS_complete == "YES" 
	)
	
	df[s2, "decision_description"] <- "Complete ITS region"
	df[s2, "Final_contig"] <- df$extraction_ITS_complete_path[s2]
	df[s2, "Final_contig_desc"] <- "complete_ITS_region"
	df[s2, "Final_outcome"] <- "PASS"
	
	### Scenario 5, 7 & 8
	s5 <- which(df$blast_its2_correct_taxonomy	== "PASS" & 
				df$blast_its1_correct_taxonomy == "PASS" &
				df$overlapping_coords == "YES" &
				df$same_contig == "YES" &
				df$extraction_ITS_complete == "NO" 
	)
	
	df[s5, "decision_description"] <- "WARNING: Overlapping coordinates in blast results"
	df[s5, "Final_contig"] <- df$blast_its2_contig_path[s5]
	df[s5, "Final_contig_desc"] <- "putative_ITS_region"
	df[s5, "Final_outcome"] <- "PASS"
	
	### Scenario 6
	s6 <- which(df$blast_its2_correct_taxonomy	== "PASS" & 
				df$blast_its1_correct_taxonomy == "PASS" &
				df$overlapping_coords == "YES" &
				df$same_contig == "YES" &
				df$extraction_ITS_complete == "YES" 
	)
	
	df[s6, "decision_description"] <- "Complete ITS region"
	df[s6, "Final_contig"] <- df$extraction_ITS_complete_path[s6]
	df[s2, "Final_contig_desc"] <- "complete_ITS_region"
	df[s6, "Final_outcome"] <- "PASS"
	
	
	### Scenario 9 & 11
	s9 <- which(df$blast_its2_correct_taxonomy	== "PASS" & 
				df$blast_its1_correct_taxonomy == "FAIL" &
				df$extraction_ITS_complete == "NO" &
				df$extraction_ITS2 == "NO"
	)
	
	df[s9, "decision_description"] <-  "Putative ITS region"
	df[s9, "Final_contig"] <- df$blast_its2_contig_path[s9]
	df[s9, "Final_contig_desc"] <- "putative_ITS_region"
	df[s9, "Final_outcome"] <- "PASS"
	
	### Scenario 10
	s10 <- which(df$blast_its2_correct_taxonomy	== "PASS" & 
				df$blast_its1_correct_taxonomy == "FAIL" &
				df$extraction_ITS_complete == "NO" &
				df$extraction_ITS2 == "YES"
	)
	
	df[s10, "decision_description"] <-  "ITS2 found only"
	df[s10, "Final_contig"] <- df$extraction_ITS2_path[s10]
	df[s10, "Final_contig_desc"] <- "ITS2"
	df[s10, "Final_outcome"] <- "PASS"
	
	### Scenario 12
	s12 <- which(df$blast_its2_correct_taxonomy	== "PASS" & 
				df$blast_its1_correct_taxonomy == "FAIL" &
				df$extraction_ITS_complete == "YES" 
	)
	
	df[s12, "decision_description"] <- "Complete ITS region"
	df[s12, "Final_contig"] <- df$df$extraction_ITS_complete_path[s12]
	df[s12, "Final_contig_desc"] <- "complete_ITS_region"
	df[s12, "Final_outcome"] <- "PASS"
	
	### Scenario 13, 15 & 16
	s13 <- which(df$blast_its2_correct_taxonomy	== "PASS" & 
				df$blast_its1_correct_taxonomy == "PASS" &
				df$same_contig == "NO" &
				df$extraction_ITS_complete == "NO" 
	)
	
	df[s13, "decision_description"] <- "WARNING: Blast results not contiguous - defaulting to ITS2 contig"
	df[s13, "Final_contig"] <- df$blast_its2_contig_path[s13]
	df[s13, "Final_contig_desc"] <- "putative_ITS_region"
	df[s13, "Final_outcome"] <- "PASS"
	
	### Scenario 14
	s14 <- which(df$blast_its2_correct_taxonomy	== "PASS" & 
				df$blast_its1_correct_taxonomy == "PASS" &
				df$same_contig == "NO" &
				df$extraction_ITS_complete == "YES" 
	)
	
	df[s14, "decision_description"] <- "Complete ITS region"
	df[s14, "Final_contig"] <- df$extraction_ITS_complete_path[s14]
	df[s14, "Final_contig_desc"] <- "complete_ITS_region"
	df[s14, "Final_outcome"] <- "PASS"
	
	
	### Scenario 17 18 & 19
	s17 <- which(df$blast_round1_correct_taxonomy == "PASS" &
				df$blast_its2_correct_taxonomy == "FAIL" &
				df$blast_its1_correct_taxonomy == "PASS" 
	)
	
	df[s17, "decision_description"] <- "Putative ITS region"
	df[s17, "Final_contig"] <- df$blast_its1_contig_path[s17]
	df[s17, "Final_contig_desc"] <- "putative_ITS_region"
	df[s17, "Final_outcome"] <- "PASS"
	
	s18 <- which(df$blast_round1_correct_taxonomy == "PASS" &
				df$blast_its2_correct_taxonomy == "FAIL" &
				df$blast_its1_correct_taxonomy == "FAIL" 
	)
	df[s18, "decision_description"] <- "Not enough support - manual curation required"
	df[s18, "Final_contig"] <- df$blast_round1_contig_path[s18]
	df[s18, "Final_outcome"] <- "FAIL"
	
	s19 <- which(df$blast_round1_correct_taxonomy == "FAIL")
	
	df[s19, "decision_description"] <- "Failed all checks"
	df[s19, "Final_contig"] <- NA
	df[s19, "Final_outcome"] <- "FAIL"


        #last scenario - failure before blast rounds
        s20 <- which(df$Contig_desc == "FAILED CONTIG" & is.na(df$Final_outcome))
        df[s20, "decision_description"] <- "Failed all checks"
        df[s20, "Final_contig"] <- NA
        df[s20, "Final_outcome"] <- "FAIL"

	df
}
