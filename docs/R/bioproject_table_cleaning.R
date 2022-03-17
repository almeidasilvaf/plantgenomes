
#----Load plant ontology--------------------------------------------------------
po <- readr::read_tsv("https://raw.githubusercontent.com/Planteome/plant-ontology/master/plant-ontology.txt",
                      skip = 1, show_col_types = FALSE)
po <- po[!startsWith(po$defn, "OBSOLETE"), ]

#----Load bioproject table------------------------------------------------------
brassicaceae_path <- paste0(tempdir(), "/brassicaceae.rda")
poaceae_path <- paste0(tempdir(), "/poaceae.rda")
fabaceae_path <- paste0(tempdir(), "/fabaceae.rda")
otherfam_path <- paste0(tempdir(), "/otherfam.rda")

# Download files in tempdir()
download.file("https://github.com/almeidasilvaf/plant_transcriptomes/blob/master/data/brassicaceae_projects.rda?raw=true", destfile = brassicaceae_path)
download.file("https://github.com/almeidasilvaf/plant_transcriptomes/blob/master/data/poaceae_projects.rda?raw=true", destfile = poaceae_path)
download.file("https://github.com/almeidasilvaf/plant_transcriptomes/blob/master/data/fabaceae_projects.rda?raw=true", destfile = fabaceae_path)
download.file("https://github.com/almeidasilvaf/plant_transcriptomes/blob/master/data/otherfam_projects.rda?raw=true", destfile = otherfam_path)

# Read files
load(brassicaceae_path)
load(fabaceae_path)
load(poaceae_path)
load(otherfam_path)

bioprojects <- bind_rows(
    brassicaceae_projects %>% dplyr::distinct(., .keep_all = TRUE),
    poaceae_projects %>% dplyr::distinct(., .keep_all = TRUE),
    fabaceae_projects %>% dplyr::distinct(., .keep_all = TRUE),
    otherfam_projects %>% dplyr::distinct(., .keep_all = TRUE),
) %>%
    mutate(BioProject = as.factor(BioProject),
           N = as.numeric(N),
           Pubmed = as.factor(Pubmed),
           Species = as.factor(Species))


#----Standardize tissue names---------------------------------------------------
# See all unique entries for the 'Tissue' variable and their frequencies
check_tissue <- function(bp = NULL) {
    tissues <- table(bp$tissue_name)
    tissue_df <- data.frame(tissues)
    tissue_df <- tissue_df[order(tissue_df$Freq, decreasing = TRUE), ]
    names(tissue_df)[1] <- "Tissue" 
    return(tissue_df)
}

# Check if tissue is included in Plant Ontology
check_exists <- function(tissue = NULL) {
    e <- tissue %in% po$name
    if(e) {
        e <- po$id[po$name == tissue]
        e <- e[!is.na(e)]
        e <- paste0(tissue, " (", e, ")")
    }
    return(e)
}


# Data frame of bioprojects containing standard Plant Ontology Terms
bp <- bioprojects %>%
    separate_rows(Tissue, sep = " \\| ", convert = TRUE) %>%
    separate(Tissue, c("tissue_name", "tissue_n"), sep = ": ") %>%
    mutate(tissue_name = str_to_title(tissue_name)) %>%
    mutate(tissue_name = str_replace_all(
        tissue_name,
        c(".*[Ss]ource_name" = "",
          " Type" = "",
          "_Type" = "",
          ".*[Mm]esocotyl.*" = "mesocotyl (PO:0020037)",
          ".*[Ss]ilk.*" = "silk (PO:0006488)",
          ".*[Oo]var[yi].*" = "plant ovary (PO:0009072)",
          ".*[Ff]iber.*" = "plant fiber cell (PO:0025407)",
          ".*[Pp]hloem.*" = "phloem (PO:0005417)",
          ".*[Mm]eiocyte.*" = "sporocyte (PO:0006204)",
          ".*[Tt]richome.*" = "trichome (PO:0000282)",
          ".*[Pp]etiole.*" = "petiole (PO:0020038)",
          ".*[Xx]ylem.*" = "xylem (PO:0005352)",
          ".*[Bb]ranch.*" = "branch (PO:0025073)",
          ".*[Ee]picotyl.*" = "epicotyl (PO:0020035)",
          ".*[Ee]xocarp.*" = "exocarp (PO:0009085)",
          ".*[Ff]oliage.*" = "",
          ".*[Hh]usk.*" = "inflorescence bract (PO:0009054)",
          ".*[Ff]oliage.*" = "leaf (PO:0025034)",
          ".*[Jj]uice Sac.*" = "juice sac tissue (PO:0006013)",
          ".*[Pp]lacenta.*" = "peduncle (PO:0009053)",
          ".*[Pp]eduncle.*" = "peduncle (PO:0009053)",
          ".*[Rr]adicle.*" = "radicle (PO:0020031)",
          ".*[Aa]ril.*" = "aril (PO:0009090)",
          ".*[Cc]ob.*" = "ear infructescence axis (PO:0025623)",
          ".*[Cc]oleoptile.*" = "seedling coleoptile (PO:0025287)",
          ".*[Ff]ilament.*" = "filament (PO:0009067)",
          ".*[Cc]ollenchyma.*" = "collenchyma (PO:0005423)",
          ".*[Hh]ypanthi.*" = "hypanthium (PO:0009065)",
          ".*[Nn]ucell.*" = "nucellus (PO:0020020)",
          ".*[Mm]ale Gametophyte.*" = "microgametophyte (PO:0025280)",
          ".*Panciles.*" = "panicle inflorescence (PO:0030123)",
          ".*[Pp]ith.*" = "pith (PO:0006109)",
          ".*[Ll]igule.*" = "ligule (PO:0020105)",
          ".*[Ss]tipule.*" = "stipule (PO:0020041)",
          ".*[Ss]uspensor.*" = "suspensor (PO:0020108)",
          ".*[Tt]esta.*" = "testa (PO:0020057)",
          ".*[Tt]hall.*" = "thallus (PO:0030027)",
          ".*[Tt]orus.*" = "pith (PO:0006109)",
          ".*[Vv]ein.*" = "vascular bundle (PO:0005020)",
          ".*[Ww]ood.*" = "stem (PO:0009047)",
          ".*[Aa]leurone.*" = "aleurone layer (PO:0005360)",
          ".*[Cc]hlorenchyma.*" = "chlorenchyma (PO:0005426)",
          ".*[Aa]rchegoni.*" = "archegonium (PO:0025126)",
          ".*[Cc]olumella.*" = "columella (PO:0025230)",
          ".*[Cc]alyx.*" = "calyx (PO:0009060)",
          ".*[Gg]ynoeci.*" = "gynoecium (PO:0009062)",
          ".*[Gg]ynophore.*" = "gynophore (PO:0006330)",
          ".*[Pp]lumul.*" = "plumule (PO:0020032)",
          ".*[Zz]ygote.*" = "plant zygote (PO:0000423)",
          ".*[Ss]uspension Cells" = "cultured plant cell (PO:0000005)",
          ".*[Ss]tolon.*" = "stolon (PO:0003024)",
          ".*[Ee]piderm.*" = "epidermis (PO:0005679)",
          ".*[Cc]apsule.*" = "capsule fruit (PO:0030091)",
          ".*[Ff]loret.*" = "tassel floret (PO:0006310)",
          ".*[Cc]hloronema.*" = "chloronema (PO:0030004)",
          ".*[Gg]lume.*" = "glume (PO:0009039)",
          ".*[Mm]icrospore.*" = "microspore (PO:0020048)",
          ".*[Bb]ract.*" = "bract (PO:0009055)",
          ".*[Cc]aryopsis.*" = "caryopsis fruit (PO:0030104)",
          ".*[Cc]orolla.*" = "corolla (PO:0009059)",
          ".*[Cc]ell Culture" = "cultured plant cell (PO:0000005)",
          ".*[Bb]err[yi].*" = "berry fruit (PO:0030108)",
          ".*[Bb]aya.*" = "berry fruit (PO:0030108)",
          ".*[Pp]arenchyma.*" = "parenchyma (PO:0005421)",
          ".*[Dd]rupe.*" = "drupe fruit (PO:0030103)",
          ".*[Cc]ytoledon.*" = "cotyledon (PO:0020030)",
          ".*[Ee]ndocarp.*" = "endocarp (PO:0009086)",
          ".*[Ee]gg Cell.*" = "plant egg cell (PO:0020094)",
          ".*[Ee]pidemis.*" = "epidermis (PO:0005679)",
          ".*[Ee]xocarp.*" = "exocarp (PO:0009085)",
          ".*[Ee]ndorsperm.*" = "endosperm (PO:0009089)",
          ".*[Mm]esophyl.*" = "mesophyll (PO:0006070)",
          ".*[Pp]neumatophore.*" = "pneumatophore (PO:0025357)",
          ".*[Pp]rotoplast.*" = "plant protoplast (PO:0000006)",
          ".*[Ss]eptum.*" = "septum (PO:0000030)",
          ".*[Ss]cutellum.*" = "scutellum (PO:0020110)",
          ".*[Ee]ndosperm.*" = "endosperm (PO:0009089)",
          ".*[Gg]rain.*" = "seed (PO:0009010)",
          ".*[Pp]ulp.*" = "endocarp (PO:0009086)",
          ".*[Tt]assel.*" = "tassel spikelet (PO:0006309)",
          ".*[Rr]hizome.*" = "rhizome (PO:0004542)",
          ".*[Ss]pike.*" = "spike inflorescence (PO:0030117)",
          ".*[Tt]hallus.*" = "thallus (PO:0030027)",
          ".*[Ll]atex.*" = "laticifer cell (PO:0025031)",
          ".*[Bb]ark.*" = "bark (PO:0004518)",
          ".*[Cc]oleoptil.*" = "coleoptile (PO:0020033)",
          ".*[Cc]oleoptyl.*" = "coleoptile (PO:0020033)",
          ".*[Gg]ametophore.*" = "gametophore (PO:0030018)",
          ".*[Ss]taminod.*" = "staminode (PO:0009077)",
          ".*[Aa]lbedo.*" = "mesocarp (PO:0009087)",
          ".*[Pp]rotonema.*" = "protonema (PO:0030003)",
          ".*[Rr]eceptacle.*" = "receptacle (PO:0009064)",
          ".*[Ss]talk.*" = "stalk (PO:0025066)",
          ".*[Tt]endril.*" = "leaf tendril (PO:0025361)",
          ".*[Hh]ydathode.*" = "hydathode (PO:0005660)",
          ".*[Ll]eaf.*" = "leaf (PO:0025034)",
          ".*[Ll]eaves.*" = "leaf (PO:0025034)",
          ".*[Ll]eave.*" = "leaf (PO:0025034)",
          ".*[Rr]osette.*" = "leaf (PO:0025034)",
          ".*[Ss]heath.*" = "leaf (PO:0025034)",
          ".*[Bb]lade.*" = "leaf (PO:0025034)",
          ".*[Ff]lesh.*" = "fruit (PO:0009001)",
          ".*[Ff]ruit.*" = "fruit (PO:0009001)",
          ".*[Pp]eel.*" = "pericarp (PO:0009084)",
          ".*[Pp]ericarp.*" = "pericarp (PO:0009084)", # fruit
          ".*[Mm]esocarp.*" = "mesocarp (PO:0009087)",
          ".*[Ss]arcocarp.*" = "mesocarp (PO:0009087)",
          ".*[Ss]eedling.*" = "seedling development stage (PO:0007131)",
          ".*[Ss]eeding.*" = "seedling development stage (PO:0007131)",
          ".*[Rr]oot Nodule.*" = "root nodule (PO:0003023)",
          ".*[Nn]odule.*" = "root nodule (PO:0003023)",
          ".*[Ff]lower.*" = "flower (PO:0009046)",
          ".*[Ff]loral.*" = "flower (PO:0009046)",
          ".*[Cc]arpel.*" = "carpel (PO:0009030)", # flower
          ".*[Aa]nther.*" = "anther (PO:0009066)", # flower
          ".*[Pp]istil.*" = "carpel (PO:0009030)", # flower
          ".*[Pp]edicel.*" = "pedicel (PO:0030112)", # flower
          ".*[Ss]tamen.*" = "stamen (PO:0009029)", # flower
          ".*[Pp]etal.*" = "petal (PO:0009032)", # flower
          ".*[Ss]epal.*" = "sepal (PO:0009031)", # flower
          ".*[Ss]tigma.*" = "stigma (PO:0009073)", # flower,
          ".*[Ss]tyle.*" = "style (PO:0009074)", # flower
          ".*[Pp]ollen.*" = "style (PO:0009074)", # flower
          ".*[Nn]ectary.*" = "nectary (PO:0009035)", # flower
          ".*[Pp]od.*" = "legume fruit (PO:0030100)",
          ".*[Ss]hoot.*" = "shoot system (PO:0009006)",
          ".*[Ss]cion.*" = "shoot system (PO:0009006)",
          ".*[Aa]erial.*" = "shoot system (PO:0009006)",
          ".*[Ee]picarp.*" = "exocarp (PO:0009085)",
          ".*[Aa]boveground.*" = "shoot system (PO:0009006)",
          ".*[Ss]am .*" = "shoot meristematic apical cell (PO:0030009)",
          ".*[Ss]am-.*" = "shoot meristematic apical cell (PO:0030009)",
          ".*[Ss]am$" = "shoot meristematic apical cell (PO:0030009)",
          ".*[Ss]am,.*" = "shoot meristematic apical cell (PO:0030009)",
          ".*[Bb]ud.*" = "bud (PO:0000055)", # Shoot
          ".*[Aa]bove Ground.*" = "shoot system (PO:0009006)",
          ".*[Aa]pical.*" = "shoot axis apex (PO:0000037)",
          ".*[Aa]pex.*" = "shoot axis apex (PO:0000037)",
          ".*[Ss]eed.+[Cc]oat.*" = "seed coat (PO:0009088)",
          ".*[Ss]eed .*" = "seed (PO:0009010)",
          ".*[Ss]eeds.*" = "seed (PO:0009010)",
          ".*[Ss]eed$" = "seed (PO:0009010)",
          ".*[Ss]eed[0-9].*" = "seed (PO:0009010)",
          ".*[Ss]eed-.*" = "seed (PO:0009010)",
          ".*[Ss]eed_.*" = "seed (PO:0009010)",
          ".*[Ii]nflorescence.*" = "inflorescence (PO:0009049)",
          ".*[Ii]nfloresence.*" = "inflorescence (PO:0009049)",
          ".*[Ii]nflorescense.*" = "inflorescence (PO:0009049)",
          ".*[Pp]anicle.*" = "inflorescence (PO:0009049)",
          ".*[Kk]ernel.*" = "ear infructescence (PO:0025597)",
          ".*[Ee]mbryogenic.*" = "embryogenic callus (PO:0006091)",
          ".*[Ee]mbryo [Aa]xis" = "plant embryo axis (PO:0019018)",
          ".*[Ee]mbryonic [Aa]xis" = "plant embryo axis (PO:0019018)",
          ".*[Ss]omatic.*[Ee]mbryo.*" = "somatic plant embryo (PO:0025302)",
          ".*[Ee]mbryo.*" = "plant embryo proper (PO:0000001)",
          ".*[Ee]mbyro*" = "plant embryo proper (PO:0000001)",
          ".*[Rr]oot.*" = "root (PO:0009005)",
          ".*[Ss]ilique.*" = "silique fruit (PO:0030106)",
          ".*[Ss]tem.*" = "stem (PO:0009047)",
          ".*[Ii]nternode.*" = "stem (PO:0009047)",
          ".*[Oo]vule.*" = "plant ovule (PO:0020003)",
          ".*[Ee]ndosperm.*" = "endosperm (PO:0009089)",
          ".*[Cc]allus.*" = "plant callus (PO:0005052)",
          ".*[Cc]alus.*" = "plant callus (PO:0005052)",
          ".*[Cc]alli.*" = "plant callus (PO:0005052)",
          ".*[Hh]ypocotyl.*" = "hypocotyl (PO:0020100)",
          ".*[Cc]otyledon.*" = "cotyledon (PO:0020030)",
          ".*[Ee]ar.*" = "ear spikelet (PO:0006320)",
          ".*[Tt]uber.*" = "tuber (PO:0025522)",
          ".*[Pp]ool.*" = "whole plant (PO:0000003)",
          ".*[Cc]ortex.*" = "cortex (PO:0005708)",
          ".*[Ww]hole Plant.*" = "whole plant (PO:0000003)",
          ".*[Ww]hole.*" = "whole plant (PO:0000003)",
          ".*[Ss]porophyte.*" = "whole plant (PO:0000003)"
        )
    )
    ) 

bp_po <- bp %>%
    filter(str_detect(tissue_name, "\\(PO")) %>%
    filter(!str_detect(tissue_n, "[A-Za-z]")) %>%
    mutate(tissue_n = as.numeric(tissue_n)) %>%
    group_by(BioProject, tissue_name, Species) %>%
    mutate(tissue_n = sum(tissue_n)) %>%
    distinct(., .keep_all = TRUE) %>%
    ungroup() %>%
    arrange(BioProject) %>%
    group_by(BioProject, Species) %>%
    summarise(tissue_count = stringr::str_c(tissue_name, ": ", tissue_n, 
                                            collapse = " | ")) %>%
    ungroup()

bp_po <- bp %>%
    select(BioProject, N, `Study title`, `Study abstract`, Pubmed, Species) %>%
    inner_join(., bp_po) %>%
    dplyr::rename(Tissue = tissue_count) %>%
    select(BioProject, N, Tissue, 
           `Study title`, `Study abstract`, Pubmed, Species) %>%
    distinct(., .keep_all = TRUE)
    

# Data frame of bioprojects containing poorly labeled tissues - set to NA
bp_not_po <- bp %>%
    filter(!str_detect(tissue_name, "\\(PO")) %>%
    mutate(tissue_name = NA) %>%
    select(-tissue_n) %>%
    distinct(., .keep_all = TRUE) %>%
    dplyr::rename(Tissue = tissue_name)

# Data frame of bioprojects that have NA in 'Tissue' 
# (they were discarded from bp_po and bp_not_po)
bp_rest <- bioprojects[is.na(bioprojects$Tissue), ]

# Combining everything
projects <- bind_rows(bp_po, bp_not_po, bp_rest)

save(
    projects,
    file = here::here("data", "projects.rda"),
    compress = "xz"
)

