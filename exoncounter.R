#the bit that i wish i'd done in the first intronexon program.
file_paths <- list.files(pattern = "\\.gtf$")

results <- data.frame(Species = character(), Exon_Count = integer())


for (file_path in file_paths) {
  gtf_lines <- readLines(file_path)
  exon_count <- sum(grepl("exon", gtf_lines, ignore.case = TRUE)) #looks for exons, ignores cases, i dont trust gtf files to not randomly have some capital letters in there
  species_name <- sub("\\.gtf$", "", file_path) # chops off the gtf bit
  results <- rbind(results, data.frame(Species = species_name, Exon_Count = exon_count)) # throws everything into dataframe as new rows
}
write.table(results, "exon_counts.csv", row.names = FALSE, col.names = c("Species", "Exon_Count"), sep = ",", quote = FALSE) #writes to file :)

# print the results because im too lazy to go and open the csv
print(results)
