#Scrip to plot a treemap of miRNA-disease association data
#   from a dataframe of microRNAs presenting SNPs

## load libraries
library("dplyr")
library("stringr")
library("vroom")
library("treemap")

## Read args from command line
args = commandArgs(trailingOnly=TRUE)

## Uncomment For debugging only
## Comment for production mode only
#args[1] <-"./test/data/sample_out.tsv" # TSV test file with miRNA data

#args[2] <- "sample" # output png file

## get the BED files
miRNAs <- args[1]

## pass to named objects
output_baseName <- args[2]

##Read files

# Read miRNA-disease association data from HMDD (the Human microRNA Disease Database)
# https://www.cuilab.cn/static/hmdd3/data/alldata.txt
HMDD_data <- vroom(file = "HMDD_alldata.txt")

## Uncomment For debugging only
#HMDD_data <- vroom(file = "https://www.cuilab.cn/static/hmdd3/data/alldata.txt")

miRNAs.df <- read.table(file = miRNAs , sep = "\t", header = T )

seed_diseases.v <- miRNAs.df %>%
    filter(type == "miRNA_seed") %>%
    pull (mir) %>%
    str_replace( "R" , "r")

seed_diseases.df <- HMDD_data %>%
    filter(mir %in% seed_diseases.v )

mature_diseases.v <- miRNAs.df %>%
    filter(type == "miRNA") %>%
    pull (mir) %>%
    str_replace( "R" , "r")

mature_diseases.df <- HMDD_data %>%
    filter(mir %in% mature_diseases.v )

primary_diseases.v <- miRNAs.df %>%
    filter(type == "miRNA_primary_transcript") %>%
    pull (mir) %>%
    str_replace( "R" , "r")

primary_diseases.df <- HMDD_data %>%
    filter(mir %in% primary_diseases.v )

#Write ouput dataframe with the miRNA and disease association data
write.table(seed_diseases.df,
            file = str_interp("${output_baseName}_seed_diseases.tsv"),
            sep = "\t",
            row.names = F, col.names = T, quote = F)

write.table(mature_diseases.df,
            file = str_interp("${output_baseName}_mature_diseases.tsv"),
            sep = "\t",
            row.names = F, col.names = T, quote = F)

write.table(primary_diseases.df,
            file = str_interp("${output_baseName}_primary_diseases.tsv"),
            sep = "\t",
            row.names = F, col.names = T, quote = F)


#Count disease association data and print top 10 more frequent with miRNA region data
count_primary.df <- primary_diseases.df %>%
    count(disease, sort = TRUE) %>%  head(10)


count_mature.df <- mature_diseases.df %>%
    count(disease, sort = TRUE) %>%  head(10)


count_seed.df <- seed_diseases.df %>%
    count(disease, sort = TRUE)

# Plot treemap with top 10 disease association data according to miRNA region data
# treemap
if(nrow(count_primary.df)>0){
    png(filename= str_interp("${output_baseName}_primary.png"),
        units="in",
        width=8,
        height=8,
        pointsize=12,
        res=72)

    treemap(count_primary.df,
            index="disease",
            vSize="n",
            type="index",
            palette = "Reds",
            title="Diseases associated with primary miRNAs that presented SNVs", #Customize your title
            fontsize.title = 12 #Change the font size of the title
    )
    dev.off()
}
# Seed region
if(nrow(count_seed.df)>0){
    png(filename= str_interp("${output_baseName}_seed.png"),
        units="in",
        width=8,
        height=8,
        pointsize=12,
        res=72)

    treemap(count_seed.df,
            index="disease",
            vSize="n",
            type="index",
            palette = "Reds",
            title="Diseases associated with miRNAs that presented SNVs in seed region", #Customize your title
            fontsize.title = 09 #Change the font size of the title
    )

    dev.off()
}

# Mature region
if(nrow(count_mature.df)>0){
    png(filename= str_interp("${output_baseName}_mature.png"),
        units="in",
        width=8,
        height=8,
        pointsize=12,
        res=72)

    treemap(count_mature.df,
            index="disease",
            vSize="n",
            type="index",
            palette = "Reds",
            title="Diseases associated with mature miRNAs that presented SNVs", #Customize your title
            fontsize.title = 09 #Change the font size of the title
    )

    dev.off()
}

