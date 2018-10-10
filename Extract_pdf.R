# This code extracts pdfs from a list of hyperlinks
# Finished, but a few other things could be added

library(rvest); library(stringr); library(pdftools)

tweets = read.csv("~/tweets.csv")

for(i in 1:5){
  tweets$expanded_urls = gsub("(.*),.*", "\\1", tweets$expanded_urls) # run this command a couple of times if many commas are present
}


links = 1:nrow(tweets)
stopwords = c(":")
missing = as.data.frame(matrix(nrow = nrow(tweets), ncol = 1)); names(missing) = "Missing.pdfs"
paper = 1


for(i in links){
  page = try(html(paste("", tweets[i,10], "", sep = "")))
  
  if(isTRUE(class(page)=="try-error")) { next } 
  else
    
  page = page %>%
    html_nodes("a") %>%       # find all links
    html_attr("href") %>%     # get the url
    str_subset("\\.pdf") # find those that end in pdf
    
  page = as.list(page)
  page = grep(paste0("http", collapse="|"), page, value = TRUE)
  
  for(j in page){
    a = try(download.file(page, "test.pdf", mode = "wb"))
    
    if(isTRUE(class(a)=="try-error")) { missing[paper,] = as.character(tweets[paper,6]) } 
    else
      
    pdf_file = file.path("C:/Users/Alfred/Downloads/Twitter/Twitter-papers", "test.pdf")
    name_file = pdf_info(pdf_file)$keys$Title
    name_file = sub(paste0(stopwords, collapse = "|"), "", name_file)
    name_file = gsub("\n", " ", name_file)
    
    if(isTRUE(sapply(strsplit(name_file, " "), length) < 1)) { 
      name_file = i
      file.rename(pdf_file, paste(name_file, ".pdf", sep = ""))
    } 
    else
      
    check_supp = pdf_text(pdf_file)[1]
    pattern = str_detect(check_supp, c("Supplementary material", "supplementary material", 
                                       "Supplementary information", "supplementary information",
                                       "SUPPLEMENTARY INFORMATION", "SUPPLEMENTARY MATERIAL"))
    
    if(isTRUE(pattern[1] | pattern[2] | pattern[3] | pattern[4] | pattern[5] | pattern[6])) { 
      name_file = pdf_info(pdf_file)$keys$Author
      file.rename(pdf_file, paste(name_file, "_Supplementary_material", ".pdf", sep = ""))
    } 
    else
      
      file.rename(pdf_file, paste(name_file, ".pdf", sep = ""))
  }
  
  paper = paper + 1
  
}
