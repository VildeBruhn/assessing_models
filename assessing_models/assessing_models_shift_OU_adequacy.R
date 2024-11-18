
#---------------------------------------
# Testing the adequacy of the OU models
#---------------------------------------

#R version 4.2.1
#evoTS version 1.0.3
#adePEM new models version
#paleoTS version 0.6.2

load('./results_paleoTS_v0.6.1/Results_OUsubsets_shiftmodels.RData')

# test adequacy

OU_adeq <- list()
OU_mov_opt_anc_adeq <- list()
OU_mov_opt_adeq <- list()
Stasis_OU_subset2_adeq <- list()
URW_OU_subset2_adeq <- list()
GRW_OU_subset2_adeq <- list()
OU_OU_subset1_adeq <- list() 
OU_OU_subset2_adeq <- list()
OU_GRW_subset1_adeq <- list()
OU_URW_subset1_adeq <- list()
OU_Stasis_subset1_adeq <- list()


if (length(OU) > 0) { 
  for (i in 1:length(OU)) {
   OU_adeq[[i]] <- fit3adequacy.OU(OU[[i]], plot = FALSE)
  }
}

if (length(OU_mov_opt_anc) > 0) { 
  for (i in 1:length(OU_mov_opt_anc)) {
   OU_mov_opt_anc_adeq[[i]] <- fit3adequacy.OU(OU_mov_opt_anc[[i]], plot = FALSE)
  }
}

if (length(OU_mov_opt) > 0) { 
  for (i in 1:length(OU_mov_opt)) {
    OU_mov_opt_adeq[[i]] <- fit3adequacy.OU(OU_mov_opt[[i]], plot = FALSE)
  }
}

if (length(Stasis_OU_subset2) > 0) { 
  for (i in 1:length(Stasis_OU_subset2)) {
   Stasis_OU_subset2_adeq[[i]] <- fit3adequacy.OU(Stasis_OU_subset2[[i]], plot = FALSE) 
  }
}

if (length(URW_OU_subset2) > 0) { 
  for (i in 1:length(URW_OU_subset2)) {
    URW_OU_subset2_adeq[[i]] <- fit3adequacy.OU(URW_OU_subset2[[i]], plot = FALSE)
  }
}

if (length(GRW_OU_subset2) > 0) { 
  for (i in 1:length(GRW_OU_subset2)) {
   GRW_OU_subset2_adeq[[i]] <- fit3adequacy.OU(GRW_OU_subset2[[i]], plot = FALSE) 
  }
}

if (length(OU_OU_subset1) > 0) { 
  for (i in 1:length(OU_OU_subset1)) {
    OU_OU_subset1_adeq[[i]] <- fit3adequacy.OU(OU_OU_subset1[[i]], plot = FALSE) 
  }
}

if (length(OU_OU_subset2) > 0) { 
  for (i in 1:length(OU_OU_subset2)) {
    OU_OU_subset2_adeq[[i]] <- fit3adequacy.OU(OU_OU_subset2[[i]], plot = FALSE)
  }
}

if (length(OU_GRW_subset1) > 0) { 
  for (i in 1:length(OU_GRW_subset1)) {
    OU_GRW_subset1_adeq[[i]] <- fit3adequacy.OU(OU_GRW_subset1[[i]], plot = FALSE) 
  }
}

if (length(OU_URW_subset1) > 0) {
  for (i in 1:length(OU_URW_subset1)) {
    OU_URW_subset1_adeq[[i]] <- fit3adequacy.OU(OU_URW_subset1[[i]], plot = FALSE)
  }
}

if (length(OU_Stasis_subset1) > 0) {
  for (i in 1:length(OU_Stasis_subset1)) {
    OU_Stasis_subset1_adeq[[i]] <- fit3adequacy.OU(OU_Stasis_subset1[[i]], plot = FALSE)
  }
}



save(OU_adeq, OU_mov_opt_anc_adeq, OU_mov_opt_adeq, Stasis_OU_subset2_adeq, 
           URW_OU_subset2_adeq, GRW_OU_subset2_adeq, OU_OU_subset1_adeq, OU_OU_subset2_adeq, 
           OU_GRW_subset1_adeq, OU_URW_subset1_adeq, OU_Stasis_subset1_adeq, file="./results_paleoTS_v0.6.1/Results_OUadeq_shiftmodels.RData")

