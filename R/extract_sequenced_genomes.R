
#----Setup----------------------------------------------------------------------
library(tidyverse)
library(rvest)
library(taxize)
library(here)

# Extract from NCBI genomes
# Search NCBI genomes by taxon
ncbi_datasets_taxon <- function(taxon = NULL) {
    
    outfile <- tempfile(fileext = ".json")
    args <- c("summary genome taxon", taxon, " > ", outfile)
    
    datasets <- system2("datasets", args = args)
    
    parsed_data <- jsonlite::fromJSON(outfile)
    unlink(outfile)
    return(parsed_data)
}

# Get BioSample location from the output of \code{ncbi_datasets_taxon()}
get_biosample_locations <- function(outlist = NULL) {
    biosamples <- outlist$assemblies$assembly$biosample
    biosamples <- biosamples[!duplicated(biosamples$accession), ]
    locations <- unlist(lapply(biosamples$attributes, function(x) {
        loc <- x$value[x$name == "geo_loc_name"]
        if(length(loc) == 0) {
            loc <- NA
        } else if(length(loc) > 1) {
            loc <- paste(loc, collapse = ", ")
        }
        return(loc)
    }))
    final_df <- data.frame(
        Biosamples = biosamples$accession,
        Location = locations
    )
    return(final_df)
}

# Get family for NCBI species IDS
get_family <- function(species_ids = NULL) {
    uspecies <- unique(species_ids)
    classification <- Reduce(rbind, lapply(uspecies, function(x) {
        fam <- taxize::classification(x, db = "ncbi")[[1]]
        family <- fam$name[fam$rank == "family"]
        order <- fam$name[fam$rank == "order"]
        class <- fam$name[fam$rank == "class"]
        phylum <- fam$name[fam$rank == "phylum"]
        
        if(length(family) == 0) { family <- NA }
        if(length(order) == 0) { order <- NA }
        if(length(class) == 0) { class <- NA }
        if(length(phylum) == 0) { phylum <- NA }
        
        fam_df <- data.frame(
            Species = x, 
            Family = family,
            Order = order, 
            Class = class,
            Phylum = phylum
        )
        return(fam_df)
    }))
    
    return(classification)
}

# Get sequencing technologies for assemblies
get_seq_tech <- function(assembly_acc = NULL, verbose = FALSE) {
    tech <- unlist(lapply(assembly_acc, function(x) {
        if(verbose) {
            message("Working on genome accession ", x)
        }
        url <- paste0("https://www.ncbi.nlm.nih.gov/assembly/", x, "/") 
        page <- rvest::read_html(url)
        name <- rvest::html_text(rvest::html_elements(page, "#summary dl > dt"))
        value <- rvest::html_text(rvest::html_elements(page, "#summary dd"))
        diff <- length(value) - length(name)
        if(diff > 0) {
            remove <- seq(3, 3 + (diff-1))
            value <- value[-c(remove)]
        }
        df <- data.frame(name, value)
        
        stech <- df$value[grepl("Sequencing technology", df$name)]
        if(length(stech) == 0) {
            stech <- NA
        }
        return(stech)
    }))
    return(tech)
}


# Include ~/Documents/Programs in RStudio's PATH
Sys.setenv(PATH = paste(Sys.getenv("PATH"), "/home/faalm/Documents/Programs", 
                        sep = ":"))

#----Get NCBI summary table-----------------------------------------------------
viridiplantae <- ncbi_datasets_taxon(taxon = "Viridiplantae")$assemblies


#----Create variables-----------------------------------------------------------
# Species names
species <- data.frame(
    Name = viridiplantae$assembly$org$sci_name,
    Species = viridiplantae$assembly$org$tax_id
)

# Taxonomy for each species
taxize_options(ncbi_sleep = 1)
tax <- get_family(viridiplantae$assembly$org$tax_id) 

# Species and family
taxonomy <- full_join(species, tax) %>%
    select(Name, Family, Order, Class, Phylum) %>%
    rename(Species = Name) %>%
    drop_na()

# Technical details
buscos <- viridiplantae$assembly$annotation_metadata$busco$complete
ngenomes <- viridiplantae$assembly$org$assembly_counts$subtree 
genome_size <- as.numeric(viridiplantae$assembly$estimated_size) / 10^6

chunk2 <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE)) 
assembly_list <- chunk2(viridiplantae$assembly$assembly_accession, 5)
tech1 <- get_seq_tech(assembly_list[[1]], verbose = TRUE)
tech2 <- get_seq_tech(assembly_list[[2]], verbose = TRUE)
tech3 <- get_seq_tech(assembly_list[[3]], verbose = TRUE)
tech4 <- get_seq_tech(assembly_list[[4]], verbose = TRUE)
tech5 <- get_seq_tech(assembly_list[[5]], verbose = TRUE)
technology <- c(tech1, tech2, tech3, tech4, tech5)

technical_df <- data.frame(
    Species = species$Name,
    BUSCOs = buscos,
    N_genomes = ngenomes,
    genome_size = round(genome_size, 2),
    technology = technology
)

# Ploidy level
ploidy <- readr::read_tsv(here::here("data", "ploidy_cvalues_db.tsv"),
                          show_col_types = FALSE) %>%
    mutate(species = paste(Genus, Species, Subspecies, sep = " ")) %>%
    select(species, `Ploidy Level (x)`) %>%
    rename(Species = species, Ploidy = `Ploidy Level (x)`) %>%
    mutate(Species = str_replace_all(Species, " NA", "")) %>%
    filter(Ploidy != "-") %>%
    mutate(Ploidy = as.numeric(Ploidy))


#----Create final data frame of genome info-------------------------------------
final_df <- merge(taxonomy, technical_df)
final_df <- merge(final_df, ploidy, all.x=TRUE)

final_df <- final_df %>%
    dplyr::distinct(., .keep_all = TRUE)

final_dfl <- split(final_df, final_df$Species)
final_dfl <- lapply(final_dfl, function(x) {
    if(nrow(x) == 1) { # only one entry per species
        df <- x
    } else {
        idx_maxbusco <- which.max(x$BUSCOs) 
        if(length(idx_maxbusco) == 0) { # none of the genomes have BUSCO info
            idx_maxsize <- which.max(x$genome_size)
            df <- x[idx_maxsize, ] 
        } else {
            df <- x[idx_maxbusco, ]
        }
    }
    return(df)
})
final_df <- Reduce(rbind, final_dfl)

# If BUSCO info is not available here and it is in the table from 
# the Nature Plants paper, get it
np <- readr::read_tsv(here("data", "nature_plants_data.txt"),
                      show_col_types = FALSE) %>%
    mutate(Species = str_c(Genus, species, sep = " "),
           NP_BUSCOs = `BUSCO % complete (combined)`) %>%
    select(Species, NP_BUSCOs) %>%
    mutate(NP_BUSCOs = str_replace_all(NP_BUSCOs, ",", ".")) %>%
    mutate(NP_BUSCOs = str_replace_all(NP_BUSCOs, "%", "")) %>%
    mutate(NP_BUSCOs = as.numeric(NP_BUSCOs)) %>%
    full_join(., final_df, by = "Species") %>%
    mutate(BUSCOs = BUSCOs * 100)

max_buscos <- pmax(np$NP_BUSCOs, np$BUSCOs, na.rm = TRUE)

final_genomes <- np %>%
    mutate(BUSCO = round(max_buscos, 2)) %>%
    select(Phylum, Class, Order, Family, Species, N_genomes, BUSCO, 
           genome_size, technology, Ploidy)

# Remove entries that only contain species info and nothing else
delete <- which(rowSums(is.na(final_genomes)) >= 9)
final_genomes <- final_genomes[-delete, ]

readr::write_tsv(
    final_genomes,
    file = here("data", "genome_info_table.tsv")
)




