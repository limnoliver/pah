# get data into properly munged format for including in the package
# ignore this on build, but push to github
library(dplyr)

dat <- read.csv('tempstorage/pahCompoundInfoGLRI.csv', na.strings = c("NA", ""), stringsAsFactors = FALSE)
dat$Parameter.code <- as.character(dat$Parameter.code)
names(dat)

dat_c <- left_join(dat, dataRetrieval::parameterCdFile, by = c('Parameter.code' = 'parameter_cd'))

# separate the compounds that have schedule = SLOH
dat_c <- filter(dat_c, !(schedule %in% 'SLOH'))

dat_c <- select(dat_c, Parameter, parameter_nm, pcode = Parameter.code, casrn,
                Abbreviation, parentAlkyl, molwt, molwt_highlow,
                benzeneRings, EPApriority16, Macdonald_SQG, VM_sourceprofiles,
                sourceProfile12, creosoteProfile11, Ingersoll_09_included,
                TEC, PEC, coc_pah_max = Coc.pah.max, coc_pah_fcv = Coc.pah.fcv)


dat_c <- mutate(dat_c,
                EPApriority16 = ifelse(is.na(EPApriority16), FALSE, TRUE),
                Macdonald_SQG = ifelse(is.na(Macdonald_SQG), FALSE, TRUE),
                VM_sourceprofiles = ifelse(is.na(VM_sourceprofiles), FALSE, TRUE),
                sourceProfile12 = ifelse(sourceProfile12 %in% 'x', TRUE, FALSE),
                creosoteProfile11 = ifelse(creosoteProfile11 %in% 'x', TRUE, FALSE),
                Ingersoll_09_included = ifelse(Ingersoll_09_included %in% 1, TRUE, FALSE))

pah_compounds <- dat_c

devtools::use_data(pah_compounds, overwrite = T)
