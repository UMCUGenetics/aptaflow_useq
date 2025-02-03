#!/usr/bin/env Rscript

library("argparse")
library("dplyr")
library("tools")
library("tidyr")
library("openxlsx")

# Parameter validation
args_parser = ArgumentParser(
  prog="SELEXrpm_top1000",
  description=""
)
args_parser$add_argument("--in-rpm-csv", "-i", help="Input file must be csv as created by SELEXrpm", required=TRUE)
args_parser$add_argument("--out-file", "-o", default="top1000.xlsx", help="Output file. Default: top1000.csv")
args_parser$add_argument("--top", "-n", type="integer", default=1000, help="Output top n sequences from every round.")
args <- args_parser$parse_args(c("-i", "rpm.csv"))


# EXCEL Workbook
wb <- createWorkbook()

# Data input
df_rpm <-  read.csv(args$in_rpm_csv, sep=";", header=TRUE)
round_names <- levels(factor(df_rpm$round))
for (round_name in round_names) {
  # Data stuff
  df_round <- df_rpm %>%
    filter(round == round_name) %>%
    arrange(desc(rpm)) %>%
    mutate(rank = row_number(desc(rpm))) %>% # rpm=round(rpm, digits=0)
    select(rank, count, rpm, seq, p5_seq_p3)
  
  df_top_1000 <- df_round %>%
    top_n(-args$top, rank) %>% # select from bottom by rank
    rename("absolute"="count", "rpm"="rpm", "random region"="seq", "full sequence"="p5_seq_p3")
  df_top_1000$"analysed as" <- " "
  
  # Fill sheet with info
  sheet = addWorksheet(wb, sheetName = round_name, gridLines = FALSE)
  writeDataTable(wb, sheet = round_name, df_top_1000, startCol = "A", startRow = 1)

}

saveWorkbook(wb, args$out_file)
