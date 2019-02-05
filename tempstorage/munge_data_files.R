# get data into properly munged format for including in the package
# ignore this on build, but push to github
library(dplyr)

#########################
# pah compound metadata
#########################
dat <- read.csv('tempstorage/pahCompoundInfoGLRI.csv', na.strings = c("NA", ""), stringsAsFactors = FALSE)
dat$Parameter.code <- as.character(dat$Parameter.code)
names(dat)

dat_c <- left_join(dat, dataRetrieval::parameterCdFile, by = c('Parameter.code' = 'parameter_cd'))

# separate the compounds that have schedule = SLOH
dat_c <- filter(dat_c, !(schedule %in% 'SLOH'))

dat_c <- select(dat_c, Parameter, parameter_nm, pcode = Parameter.code, casrn,
                Abbreviation, molwt, molwt_highlow,
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

dat_c$coc_pah_fcv[dat_c$parameter_nm %in% 'c1-alkylated fluorene'] <- 611
dat_c$coc_pah_fcv[dat_c$parameter_nm %in% 'c2-alkylated fluorene'] <- 686
dat_c$coc_pah_fcv[dat_c$parameter_nm %in% 'c3-alkylated fluorene'] <- 769

pah_compounds <- filter(dat_c, !is.na(Parameter))

devtools::use_data(pah_compounds, overwrite = T)


#########################
# pah source profiles
#########################

profiles <- read.csv('tempstorage/PAH source profiles_BaldwinEtAl2016_correctedUMO1_creosote.csv')

profiles_clean <- rename(profiles, pcode = param_cd,
                         PPLT = Power.plant.emissions,
                         RESI = Residential.heating,
                         CCB1 = Coal.average,
                         COKE = Coke.oven.emissions,
                         DVEM = Diesel.vehicle,
                         GVEM = Gasoline.vehicle,
                         TUN1 = Traffic.tunnel.air,
                         VAVG = Vehicle.traffic.avg,
                         UMO1 = Used.motor.oil.1,
                         UMO2 = Used.motor.oil.2,
                         PIN1 = Pine.combustion.1,
                         PIN2 = Pine.combustion.2,
                         OAKS = Oak.combustion,
                         FOC1 = Fuel.oil.combustion,
                         TIRE = Tire.particles,
                         ASP2 = Asphalt,
                         CTR0 = CT.T0,
                         CTR45 = CT.T45,
                         CTR376 = CT.T376,
                         CTD6 = CT.dust.6,
                         CTD7 = CT.dust.7,
                         CRE2 = Creosote.railway.ties,
                         CRE4 = Creosote.product) %>%
  mutate(pcode = gsub('P', '', pcode))

compound_info <- pah::pah_compounds
source_profiles <- left_join(profiles_clean, select(compound_info, pcode, molwt, casrn))

# test if each profile sums to 1
test <- colSums(profiles_clean[, 4:26])

devtools::use_data(source_profiles, overwrite = T)

# modify sources dataframe to add CTR0, CTR45, CTR376
sources <- pah::sources %>%
  add_row(source_subcategory = 'coal tar', source_name = c('Coal tar sealant, day 0', 'Coal tar sealant, day 45', 'Coal tar sealant, day 376'),
          source_short_name = c('Coal tar seal. prod. 0d', 'Coal tar seal. prod. 45d', 'Coal tar seal. prod. 376d'),
          source_short_no_ref = c('Coal tar seal. prod. 0d', 'Coal tar seal. prod. 45d', 'Coal tar seal. prod. 376d'),
          source_abbrev = c('CTR0', 'CTR45', 'CTR376'))


devtools::use_data(sources, overwrite = T)
