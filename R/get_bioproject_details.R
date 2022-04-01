library(bears)
library(collapse)
options(collapse_mask = "manip") 

download.file(
    "https://github.com/almeidasilvaf/plantgenomes/blob/master/data/projects.rda?raw=true",
    destfile = file.path(tempdir(), "projects.rda")
)
load(file.path(tempdir(), "projects.rda"))

bioproject_ids <- unique(as.character(projects$BioProject))

# Get bioproject info
get_bp_info <- function(ids = NULL) {
    detailed_bp <- Reduce(rbind, lapply(ids, function(x) {
        term <- paste0(x, "[GPRJ]")
        detail_df <- bears::create_sample_info(term, retmax = 1)
        detail_df <- detail_df[, c("BioProject", "Instrument", 
                                   "Cultivar", "Date", "Origin")]
        return(detail_df)
    }))
    return(detailed_bp)
}

# Repeat process in groups of 1000 projects to avoid losing everything for
# connection issues
chunk2 <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE)) 

ids <- chunk2(bioproject_ids, 10)
bp1 <- get_bp_info(ids[[1]])
bp2 <- get_bp_info(ids[[2]])
bp3 <- get_bp_info(ids[[3]])
bp4 <- get_bp_info(ids[[4]])
bp5 <- get_bp_info(ids[[5]])
bp6 <- get_bp_info(ids[[6]])
bp7 <- get_bp_info(ids[[7]])
bp8 <- get_bp_info(ids[[8]])
bp9 <- get_bp_info(ids[[9]])
bp10 <- get_bp_info(ids[[10]])

bioproject_info <- rbind(
    bp1, bp2, bp3, bp4, bp5, bp6, bp7, bp8, bp9, bp10
)

readr::write_tsv(
    bioproject_info,
    file = here::here("data", "bioproject_info.tsv")
)


