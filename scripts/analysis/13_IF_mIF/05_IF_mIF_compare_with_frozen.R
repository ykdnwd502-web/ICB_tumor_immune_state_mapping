############################################################
## 05_IF_mIF_compare_with_frozen_v1.5_FINAL.R
##
## IF/mIF frozen reproduction audit
##
## FINAL RELEASE VERSION
##
## Validation layers:
##
## Step01:
##   Step14B CosMx cell-state scores
##   key = .merge_id
##
## Step02:
##   Step14A IF/state merged table
##   key = .key_raw
##
## Step03:
##   Step14A correlation summary
##   numerical identity
##
############################################################


############################################################
## Package
############################################################

suppressPackageStartupMessages({
  
  library(dplyr)
  library(readr)
  
})


############################################################
## Project configuration
############################################################


PROJECT_DIR <- Sys.getenv(
  "ICB_PROJECT_DIR"
)


if(!nzchar(PROJECT_DIR)){
  
  PROJECT_DIR <-
    "D:/ICB_resistance_project"
  
}


PROJECT_DIR <-
  normalizePath(
    PROJECT_DIR,
    winslash="/",
    mustWork=TRUE
  )



FROZEN_ROOT <- Sys.getenv(
  "ICB_FROZEN_ROOT"
)


if(!nzchar(FROZEN_ROOT)){
  
  FROZEN_ROOT <-
    file.path(
      PROJECT_DIR,
      "results",
      "tables",
      "IF_mIF_validation",
      "frozen_reference"
    )
  
}



CURRENT_DIR <-
  file.path(
    PROJECT_DIR,
    "results",
    "tables",
    "IF_mIF_validation"
  )



AUDIT_DIR <-
  file.path(
    PROJECT_DIR,
    "results",
    "audit",
    "13_IF_mIF_sequential_reproduction_check"
  )



dir.create(
  AUDIT_DIR,
  recursive=TRUE,
  showWarnings=FALSE
)



cat("\n============================================================\n")
cat("13_IF_mIF FROZEN REPRODUCTION AUDIT v1.5 FINAL\n")
cat("============================================================\n")

cat(
  "Project:",
  PROJECT_DIR,
  "\n"
)

cat(
  "Frozen root:",
  FROZEN_ROOT,
  "\n"
)



############################################################
## Input files
############################################################


## Step14B

FROZEN_STATE_FILE <-
  file.path(
    FROZEN_ROOT,
    "Step14B_CosMx_cell_state_scores_for_IF_mIF_COMPOSITE_ID_CLEAN.csv"
  )


CURRENT_STATE_FILE <-
  file.path(
    CURRENT_DIR,
    "IF_mIF_CosMx_cell_state_scores.csv"
  )



## Step14A merged table

FROZEN_MERGE_FILE <-
  file.path(
    FROZEN_ROOT,
    "Step14A_v8_IF_mIF_state_protein_merged_table_CLEAN.csv"
  )


CURRENT_MERGE_FILE <-
  file.path(
    CURRENT_DIR,
    "IF_mIF_state_protein_merged_table.csv"
  )



## Step14A correlation

FROZEN_CORR_FILE <-
  file.path(
    FROZEN_ROOT,
    "Step14A_v8_IF_mIF_protein_state_spearman_CLEAN.csv"
  )


CURRENT_CORR_FILE <-
  file.path(
    CURRENT_DIR,
    "IF_mIF_protein_state_spearman.csv"
  )



required_files <- c(
  
  FROZEN_STATE_FILE,
  CURRENT_STATE_FILE,
  
  FROZEN_MERGE_FILE,
  CURRENT_MERGE_FILE,
  
  FROZEN_CORR_FILE,
  CURRENT_CORR_FILE
  
)



missing_files <-
  required_files[
    !file.exists(required_files)
  ]



if(length(missing_files)>0){
  
  stop(
    paste(
      "Missing files:",
      paste(
        missing_files,
        collapse="\n"
      )
    ),
    call.=FALSE
  )
  
}



############################################################
## Helper function
############################################################


numeric_maxdiff <- function(
    x,
    y,
    exclude_cols=NULL
){
  
  
  common_cols <-
    intersect(
      colnames(x),
      colnames(y)
    )
  
  
  common_cols <-
    setdiff(
      common_cols,
      exclude_cols
    )
  
  
  
  numeric_cols <-
    common_cols[
      sapply(
        x[common_cols],
        is.numeric
      )
    ]
  
  
  
  if(length(numeric_cols)==0){
    
    return(0)
    
  }
  
  
  
  diff_values <-
    sapply(
      numeric_cols,
      function(col){
        
        
        xx <- x[[col]]
        
        yy <- y[[col]]
        
        
        idx <-
          complete.cases(
            xx,
            yy
          )
        
        
        if(sum(idx)==0){
          
          return(0)
          
        }
        
        
        
        max(
          abs(
            xx[idx]-yy[idx]
          )
        )
        
        
      }
    )
  
  
  max(diff_values)
  
}



############################################################
## STEP 1
## Step14B CosMx state scores
############################################################


cat("\n============================================================\n")
cat("STEP 1: Step14B CosMx cell-state scores\n")
cat("============================================================\n")



STATE_KEY <- ".merge_id"



frozen_state <-
  read.csv(
    FROZEN_STATE_FILE,
    check.names=FALSE
  )


current_state <-
  read.csv(
    CURRENT_STATE_FILE,
    check.names=FALSE
  )



stopifnot(
  
  STATE_KEY %in%
    colnames(frozen_state),
  
  STATE_KEY %in%
    colnames(current_state)
  
)



state_ids <-
  intersect(
    frozen_state[[STATE_KEY]],
    current_state[[STATE_KEY]]
  )



fs <-
  frozen_state[
    frozen_state[[STATE_KEY]] %in% state_ids,
  ]


cs <-
  current_state[
    current_state[[STATE_KEY]] %in% state_ids,
  ]



fs <-
  fs[
    order(fs[[STATE_KEY]]),
  ]


cs <-
  cs[
    order(cs[[STATE_KEY]]),
  ]



state_diff <-
  numeric_maxdiff(
    fs,
    cs,
    exclude_cols=c(
      "fov",
      "cell_ID",
      ".merge_id"
    )
  )



STEP1_PASS <-
  length(state_ids)>0 &
  state_diff <=1e-12



cat(
  "Frozen rows:",
  nrow(frozen_state),
  "\n"
)

cat(
  "Current rows:",
  nrow(current_state),
  "\n"
)

cat(
  "Shared IDs:",
  length(state_ids),
  "\n"
)

cat(
  "Maximum difference:",
  signif(state_diff,6),
  "\n"
)





############################################################
## STEP 2
## Step14A merged table
############################################################


cat("\n============================================================\n")
cat("STEP 2: Step14A IF/state merged table\n")
cat("============================================================\n")



MERGE_KEY <- ".key_raw"



frozen_merge <-
  read.csv(
    FROZEN_MERGE_FILE,
    check.names=FALSE
  )


current_merge <-
  read.csv(
    CURRENT_MERGE_FILE,
    check.names=FALSE
  )



stopifnot(
  
  MERGE_KEY %in%
    colnames(frozen_merge),
  
  MERGE_KEY %in%
    colnames(current_merge)
  
)



merge_ids <-
  intersect(
    frozen_merge[[MERGE_KEY]],
    current_merge[[MERGE_KEY]]
  )



fm <-
  frozen_merge[
    frozen_merge[[MERGE_KEY]] %in% merge_ids,
  ]



cm <-
  current_merge[
    current_merge[[MERGE_KEY]] %in% merge_ids,
  ]



fm <-
  fm[
    order(fm[[MERGE_KEY]]),
  ]


cm <-
  cm[
    order(cm[[MERGE_KEY]]),
  ]



merge_diff <-
  numeric_maxdiff(
    fm,
    cm,
    exclude_cols=MERGE_KEY
  )



STEP2_PASS <-
  length(merge_ids)>0 &
  merge_diff <=1e-12



cat(
  "Frozen rows:",
  nrow(frozen_merge),
  "\n"
)

cat(
  "Current rows:",
  nrow(current_merge),
  "\n"
)

cat(
  "Shared keys:",
  length(merge_ids),
  "\n"
)

cat(
  "Maximum difference:",
  signif(merge_diff,6),
  "\n"
)





############################################################
## STEP 3
## Step14A correlation
############################################################


cat("\n============================================================\n")
cat("STEP 3: Step14A correlation summary\n")
cat("============================================================\n")



frozen_corr <-
  read.csv(
    FROZEN_CORR_FILE,
    check.names=FALSE
  )


current_corr <-
  read.csv(
    CURRENT_CORR_FILE,
    check.names=FALSE
  )



corr_diff <-
  numeric_maxdiff(
    frozen_corr,
    current_corr
  )



STEP3_PASS <-
  corr_diff <=1e-12



cat(
  "Maximum difference:",
  signif(corr_diff,6),
  "\n"
)




############################################################
## Audit output
############################################################


audit_table <-
  data.frame(
    
    group=c(
      "Step01_state_scoring",
      "Step02_merge_table",
      "Step03_correlation"
    ),
    
    
    item=c(
      
      "Step14B CosMx state scores reproduce frozen output",
      
      "Step14A IF/state merged table reproduce frozen output",
      
      "Step14A correlation reproduce frozen output"
      
    ),
    
    
    
    observed=c(
      
      paste0(
        "frozen_rows=",
        nrow(frozen_state),
        ";current_rows=",
        nrow(current_state),
        ";shared_ids=",
        length(state_ids),
        ";maxdiff=",
        signif(state_diff,6)
      ),
      
      
      paste0(
        "frozen_rows=",
        nrow(frozen_merge),
        ";current_rows=",
        nrow(current_merge),
        ";shared_keys=",
        length(merge_ids),
        ";maxdiff=",
        signif(merge_diff,6)
      ),
      
      
      paste0(
        "maxdiff=",
        signif(corr_diff,6)
      )
      
    ),
    
    
    
    expected=c(
      
      "shared .merge_id; maxdiff<=1e-12",
      
      "shared .key_raw; maxdiff<=1e-12",
      
      "maxdiff<=1e-12"
      
    ),
    
    
    
    pass=c(
      
      STEP1_PASS,
      
      STEP2_PASS,
      
      STEP3_PASS
      
    ),
    
    
    
    stringsAsFactors=FALSE
    
  )



AUDIT_FILE <-
  file.path(
    AUDIT_DIR,
    "13_IF_mIF_frozen_comparison_v1.5_FINAL.csv"
  )



write.csv(
  audit_table,
  AUDIT_FILE,
  row.names=FALSE
)




############################################################
## Final summary
############################################################


cat("\n============================================================\n")

cat(
  "Step01 state scoring:",
  ifelse(STEP1_PASS,"PASS","FAIL"),
  "\n"
)


cat(
  "Step02 merge:",
  ifelse(STEP2_PASS,"PASS","FAIL"),
  "\n"
)


cat(
  "Step03 correlation:",
  ifelse(STEP3_PASS,"PASS","FAIL"),
  "\n"
)



blocking <-
  sum(
    !audit_table$pass
  )


cat(
  "Blocking failures:",
  blocking,
  "\n"
)



if(blocking==0){
  
  cat(
    "Overall: PASS\n"
  )
  
}else{
  
  cat(
    "Overall: FAIL\n"
  )
  
}



cat(
  "Audit file:",
  AUDIT_FILE,
  "\n"
)


cat("============================================================\n")