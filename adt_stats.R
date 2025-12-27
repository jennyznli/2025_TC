# FIG 1 

# ========================
# LOAD DATA
# ========================
source("config.R")

ss <- read_excel(file.path("ss", ADT_META))

ss$Lymph_Node <- as.factor(ss$Lymph_Node)
ss$Chronological_Age <- as.numeric(ss$Chronological_Age)
ss$M <- as.factor(ss$M)
ss$N <- as.factor(ss$N)
ss$Probe_Success_Rate <- as.numeric(ss$Probe_Success_Rate)
ss$IC_EpiDISH <- as.numeric(ss$IC_EpiDISH)

betas <- readRDS(file.path("data", ADT_BETAS))[, ss$IDAT]
ss$Global_Means <- colMeans(betas)

# ========================
# LOAD DATA
# ========================
primary <- ss %>% filter(Lymph_Node == "F")

primary %>% filter(Lymph_Node == "F") %>% select(Sex) %>% table()
primary %>% filter(Lymph_Node == "F") %>% select(Clinical_Invasiveness) %>% table()

primary  %>% filter(Lymph_Node == "F") %>%
  summarize(
    mean = mean(Chronological_Age),
    sd = sd(Chronological_Age, na.rm = TRUE),
    min = min(Chronological_Age),
    max = max(Chronological_Age)
  )

primary %>%
  dplyr::summarize(
    n = sum(!is.na(Horvath_MethylClock)),
    median = median(Horvath_MethylClock, na.rm = TRUE),
    q1 = quantile(Horvath_MethylClock, 0.25, na.rm = TRUE),
    q3 = quantile(Horvath_MethylClock, 0.75, na.rm = TRUE),
    min = min(Horvath_MethylClock, na.rm = TRUE),
    max = max(Horvath_MethylClock, na.rm = TRUE),
    .groups = "drop"
  )


ss %>% filter(Lymph_Node == "F") %>% 
  summarize(
    mean = mean(Horvath_MethylClock),
    sd = sd(Horvath_MethylClock, na.rm = TRUE),
    min = min(Horvath_MethylClock), 
    max = max(Horvath_MethylClock)
  )

primary %>%
  dplyr::summarize(
    n = sum(!is.na(Chronological_Age)),
    median = median(Chronological_Age, na.rm = TRUE),
    q1 = quantile(Chronological_Age, 0.25, na.rm = TRUE),
    q3 = quantile(Chronological_Age, 0.75, na.rm = TRUE),
    min = min(Chronological_Age, na.rm = TRUE),
    max = max(Chronological_Age, na.rm = TRUE),
    .groups = "drop"
  )


