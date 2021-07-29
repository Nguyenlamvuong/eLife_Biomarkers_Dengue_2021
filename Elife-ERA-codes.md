Combination of inflammatory and vascular markers in the febrile phase of
dengue is associated with more severe outcomes
================
Nguyen Lam Vuong
29-Jul-2021

``` r
library(tidyverse)
library(gtsummary) # need to install package 'flextable'
library(rms)
library(MuMIn) # for best subset selection
library(facetscales) # need to install package 'facetscales' from devtools::install_github("zeehio/facetscales")
source("Elife ERA functions.R") # to include my functions
options(gtsummary.tbl_summary.percent_fun = function(x) sprintf(x * 100, fmt='%1.0f')) # to report percentages without decimal
options(knitr.kable.NA = '') # to set NA to '' in kable results
theme_set(theme_bw()) # I love black & white theme
options(na.action = "na.fail") # for 'dredge' function [MuMIn]
```

``` r
# Full data
dat0 <- read_csv("Dengue_Biomarkers_data_27Jul2021.csv") %>%
  mutate(group2 = factor(sev.or.inte, levels=c(0,1), labels=c("Uncomplicated dengue","Severe/moderate dengue")),
         Country = as.factor(Country),
         Serotype = as.factor(Serotype),
         Serology = factor(Serology, levels = c("Probable primary", "Probable secondary", "Inconclusive")),
         WHO2009 = factor(WHO2009, levels = c("Mild dengue", "Dengue with warning signs", "Severe dengue", "Unknown")))

# Data at enrollment (for Table 1 & Appendix 4-table 1)
dat <- dat0 %>% filter(Time == "Enrolment")

# Data at enrollment for models
## with inverse probability weights (IPW) for inclusion probability by countries (for the analysis of secondary outcomes)
## transform the biomarker's levels to log-2 and viremia to log-10
## set all biomarker's levels under the limit of detection (u...=1) to the limit of detection
dat1 <- dat %>%
  mutate(age15 = factor(age15, levels=c("No","Yes"), labels=c("Under 15","15 and above")),
         Serotype = ifelse(Serotype=="Unknown", NA, as.character(Serotype)),
         Serotype = ifelse(Serotype=="DENV-1", 1, 2),
         Serotype2 = Serotype, # for getting results for Appendix5-tables 1, 2
         Serotype = factor(Serotype, levels=c(1,2), labels=c("DENV-1","Others")), # set Serotype to DENV-1 and others
         ipw = ifelse(sev.or.inte==1, 1,
                      ifelse(Country=="Vietnam", (1505-204)/436,
                             ifelse(Country=="Malaysia", (259-29)/58,
                                    ifelse(Country=="El Salvador", (306-18)/23,
                                           ifelse(Country=="Cambodia", (302-30)/39, NA))))),
         ipwsd = ipw * 837/2372,
         VCAM = ifelse(uVCAM==1, log2(0.028), log2(VCAM1)),
         SDC = log2(SDC1),
         Ang = ifelse(uAng==1, log2(17.1), log2(Ang2)),
         IL8 = ifelse(uIL8==1, log2(1.8), log2(IL8)),
         IP10 = ifelse(uIP10==1, log2(1.18), log2(IP10)),
         IL1RA = ifelse(uIL1RA==1, log2(18), log2(IL1RA)),
         CD163 = log2(CD163),
         TREM = ifelse(uTREM==1, log2(10.65), log2(TREM1)),
         Fer = log2(Fer),
         CRP = log2(CRP), 
         Vir = log10(Viremia)) %>%
  select(Code, Age, age15, Serotype, Serotype2, sev.or.inte, sev.only, sev.or.ws, hospital, VCAM, SDC, Ang, IL8, IP10, IL1RA, CD163, TREM, Fer, CRP, Vir, uVCAM, uAng, uIL8, uIP10, uCD163, uTREM, uVir, ipw, ipwsd)

dat1c <- dat1 %>% filter(age15 == "Under 15") ## Data for children
dat1a <- dat1 %>% filter(age15 == "15 and above") ## Data for adults

# Set references for calculating the ORs and 95% CIs
ref0 <- c(median(dat1$VCAM), median(dat1$SDC), median(dat1$Ang), median(dat1$IL8), median(dat1$IP10), median(dat1$IL1RA), median(dat1$CD163), median(dat1$TREM), median(dat1$Fer), median(dat1$CRP))

# Create data for plotting biomarkers (for Figure 2 & Appendix 4-figure 1)
## Enrollment
tmp1 <- dat0 %>%
  filter(Time == "Enrolment") %>%
  select(Code, group2, Day, Daygr, VCAM1, SDC1, Ang2, IL8, IP10, IL1RA, CD163, TREM1, Fer, CRP) %>%
  gather(., "Biomarker", "Result", 5:14)

## Follow up
tmp2 <- dat0 %>%
  filter(Time == "Follow up") %>%
  select(Code, group2, Day, Daygr, VCAM1, SDC1, Ang2, IL8, IP10, IL1RA, CD163, TREM1, Fer, CRP) %>%
  gather(., "Biomarker", "Result", 5:14)

## Merge long data for biomarkers
dat_plot <- bind_rows(tmp1, tmp2) %>%
  mutate(Daygr = factor(Daygr, levels = c("Day 1", "Day 2", "Day 3", "Day 10-20", "Day >20"),
                        labels = c("1", "2", "3", "10-20", ">20")),
         Biomarker = factor(Biomarker, levels=c("VCAM1", "SDC1", "Ang2", "IL8", "IP10", "IL1RA", "CD163", "TREM1", "Fer", "CRP"), 
                            labels=c("VCAM-1 (ng/ml)", "SDC-1 (pg/ml)", "Ang-2 (pg/ml)", "IL-8 (pg/ml)", "IP-10 (pg/ml)", 
                                     "IL-1RA (pg/ml)", "sCD163 (ng/ml)", "sTREM-1 (pg/ml)", "Ferritin (ng/ml)", "CRP (mg/l)")))
```

########################################################################################### 

#### Table 1. Summary of clinical data by primary outcome.

``` r
# All patients
t1 <- dat %>%
  select(group2, Country, Age, Sex, Day, Serotype, Serology, Obesity, Diabetes, WHO2009, hospital) %>%
  tbl_summary(by = group2,
              statistic = list(all_continuous() ~ "{median} ({p25}, {p75})", all_categorical() ~ "{n} ({p})"),
              value = list(Sex ~ "Male"),
              digits = list(all_continuous() ~ c(0,0)),
              label = list(Age ~ "Age (years)", Sex ~ "Gender male", Day ~ "Illness day at enrolment",
                           Serology ~ "Immune status", WHO2009 ~ "WHO 2009 classification", hospital ~ "Hospitalization")) %>%
  add_stat_label(label = c(all_categorical() ~ "n (%)", Age ~ "median (1st, 3rd quartiles)")) %>%
  modify_header(label = "", stat_by = "**{level} (N={n})**")

# Children (<15 years of age)
t2 <- dat %>%
  filter(age15=="No") %>%
  select(group2, Country, Age, Sex, Day, Serotype, Serology, Obesity, Diabetes, WHO2009, hospital) %>%
  tbl_summary(by = group2,
              statistic = list(all_continuous() ~ "{median} ({p25}, {p75})", all_categorical() ~ "{n} ({p})"),
              value = list(Sex ~ "Male"),
              digits = list(all_continuous() ~ c(0,0)),
              label = list(Age ~ "Age (years)", Sex ~ "Gender male", Day ~ "Illness day at enrolment",
                           Serology ~ "Immune status", WHO2009 ~ "WHO 2009 classification", hospital ~ "Hospitalization")) %>%
  add_stat_label(label = c(all_categorical() ~ "n (%)", Age ~ "median (1st, 3rd quartiles)")) %>%
  modify_header(label = "", stat_by = "**{level} (N={n})**")

# Adults (>=15 years of age)
t3 <- dat %>%
  filter(age15=="Yes") %>%
  select(group2, Country, Age, Sex, Day, Serotype, Serology, Obesity, Diabetes, WHO2009, hospital) %>%
  tbl_summary(by = group2,
              statistic = list(all_continuous() ~ "{median} ({p25}, {p75})", all_categorical() ~ "{n} ({p})"),
              value = list(Sex ~ "Male"),
              digits = list(all_continuous() ~ c(0,0)),
              label = list(Age ~ "Age (years)", Sex ~ "Gender male", Day ~ "Illness day at enrolment",
                           Serology ~ "Immune status", WHO2009 ~ "WHO 2009 classification", hospital ~ "Hospitalization")) %>%
  add_stat_label(label = c(all_categorical() ~ "n (%)", Age ~ "median (1st, 3rd quartiles)")) %>%
  modify_header(label = "", stat_by = "**{level} (N={n})**")

# Merge tables
tbl_merge(tbls = list(t1, t2, t3),
          tab_spanner = c("**All patients**", "**Children**", "**Adults**"))
```

<div id="mfiakbcmth" style="overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>html {
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Helvetica Neue', 'Fira Sans', 'Droid Sans', Arial, sans-serif;
}

#mfiakbcmth .gt_table {
  display: table;
  border-collapse: collapse;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 16px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: auto;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}

#mfiakbcmth .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#mfiakbcmth .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#mfiakbcmth .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 0;
  padding-bottom: 4px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#mfiakbcmth .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#mfiakbcmth .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#mfiakbcmth .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}

#mfiakbcmth .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}

#mfiakbcmth .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#mfiakbcmth .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#mfiakbcmth .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}

#mfiakbcmth .gt_group_heading {
  padding: 8px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
}

#mfiakbcmth .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}

#mfiakbcmth .gt_from_md > :first-child {
  margin-top: 0;
}

#mfiakbcmth .gt_from_md > :last-child {
  margin-bottom: 0;
}

#mfiakbcmth .gt_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}

#mfiakbcmth .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 12px;
}

#mfiakbcmth .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#mfiakbcmth .gt_first_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
}

#mfiakbcmth .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#mfiakbcmth .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#mfiakbcmth .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#mfiakbcmth .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#mfiakbcmth .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#mfiakbcmth .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding: 4px;
}

#mfiakbcmth .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#mfiakbcmth .gt_sourcenote {
  font-size: 90%;
  padding: 4px;
}

#mfiakbcmth .gt_left {
  text-align: left;
}

#mfiakbcmth .gt_center {
  text-align: center;
}

#mfiakbcmth .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#mfiakbcmth .gt_font_normal {
  font-weight: normal;
}

#mfiakbcmth .gt_font_bold {
  font-weight: bold;
}

#mfiakbcmth .gt_font_italic {
  font-style: italic;
}

#mfiakbcmth .gt_super {
  font-size: 65%;
}

#mfiakbcmth .gt_footnote_marks {
  font-style: italic;
  font-weight: normal;
  font-size: 65%;
}
</style>
<table class="gt_table">
  
  <thead class="gt_col_headings">
    <tr>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="2" colspan="1"></th>
      <th class="gt_center gt_columns_top_border gt_column_spanner_outer" rowspan="1" colspan="2">
        <span class="gt_column_spanner"><strong>All patients</strong></span>
      </th>
      <th class="gt_center gt_columns_top_border gt_column_spanner_outer" rowspan="1" colspan="2">
        <span class="gt_column_spanner"><strong>Children</strong></span>
      </th>
      <th class="gt_center gt_columns_top_border gt_column_spanner_outer" rowspan="1" colspan="2">
        <span class="gt_column_spanner"><strong>Adults</strong></span>
      </th>
    </tr>
    <tr>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1"><strong>Uncomplicated dengue (N=556)</strong></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1"><strong>Severe/moderate dengue (N=281)</strong></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1"><strong>Uncomplicated dengue (N=337)</strong></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1"><strong>Severe/moderate dengue (N=127)</strong></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1"><strong>Uncomplicated dengue (N=219)</strong></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1"><strong>Severe/moderate dengue (N=154)</strong></th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr><td class="gt_row gt_left">Country, n (%)</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center"></td></tr>
    <tr><td class="gt_row gt_left" style="text-align: left; text-indent: 10px;">Cambodia</td>
<td class="gt_row gt_center">39 (7)</td>
<td class="gt_row gt_center">30 (11)</td>
<td class="gt_row gt_center">37 (11)</td>
<td class="gt_row gt_center">29 (23)</td>
<td class="gt_row gt_center">2 (1)</td>
<td class="gt_row gt_center">1 (1)</td></tr>
    <tr><td class="gt_row gt_left" style="text-align: left; text-indent: 10px;">El Salvador</td>
<td class="gt_row gt_center">23 (4)</td>
<td class="gt_row gt_center">18 (6)</td>
<td class="gt_row gt_center">23 (7)</td>
<td class="gt_row gt_center">18 (14)</td>
<td class="gt_row gt_center">0 (0)</td>
<td class="gt_row gt_center">0 (0)</td></tr>
    <tr><td class="gt_row gt_left" style="text-align: left; text-indent: 10px;">Malaysia</td>
<td class="gt_row gt_center">58 (10)</td>
<td class="gt_row gt_center">29 (10)</td>
<td class="gt_row gt_center">3 (1)</td>
<td class="gt_row gt_center">1 (1)</td>
<td class="gt_row gt_center">55 (25)</td>
<td class="gt_row gt_center">28 (18)</td></tr>
    <tr><td class="gt_row gt_left" style="text-align: left; text-indent: 10px;">Vietnam</td>
<td class="gt_row gt_center">436 (78)</td>
<td class="gt_row gt_center">204 (73)</td>
<td class="gt_row gt_center">274 (81)</td>
<td class="gt_row gt_center">79 (62)</td>
<td class="gt_row gt_center">162 (74)</td>
<td class="gt_row gt_center">125 (81)</td></tr>
    <tr><td class="gt_row gt_left">Age (years), median (1st, 3rd quartiles)</td>
<td class="gt_row gt_center">12 (9, 22)</td>
<td class="gt_row gt_center">16 (10, 24)</td>
<td class="gt_row gt_center">10 (8, 12)</td>
<td class="gt_row gt_center">10 (7, 12)</td>
<td class="gt_row gt_center">26 (20, 34)</td>
<td class="gt_row gt_center">22 (18, 30)</td></tr>
    <tr><td class="gt_row gt_left">Gender male, n (%)</td>
<td class="gt_row gt_center">299 (54)</td>
<td class="gt_row gt_center">170 (60)</td>
<td class="gt_row gt_center">173 (51)</td>
<td class="gt_row gt_center">80 (63)</td>
<td class="gt_row gt_center">126 (58)</td>
<td class="gt_row gt_center">90 (58)</td></tr>
    <tr><td class="gt_row gt_left">Illness day at enrolment, n (%)</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center"></td></tr>
    <tr><td class="gt_row gt_left" style="text-align: left; text-indent: 10px;">1</td>
<td class="gt_row gt_center">91 (16)</td>
<td class="gt_row gt_center">49 (17)</td>
<td class="gt_row gt_center">57 (17)</td>
<td class="gt_row gt_center">25 (20)</td>
<td class="gt_row gt_center">34 (16)</td>
<td class="gt_row gt_center">24 (16)</td></tr>
    <tr><td class="gt_row gt_left" style="text-align: left; text-indent: 10px;">2</td>
<td class="gt_row gt_center">260 (47)</td>
<td class="gt_row gt_center">130 (46)</td>
<td class="gt_row gt_center">156 (46)</td>
<td class="gt_row gt_center">52 (41)</td>
<td class="gt_row gt_center">104 (47)</td>
<td class="gt_row gt_center">78 (51)</td></tr>
    <tr><td class="gt_row gt_left" style="text-align: left; text-indent: 10px;">3</td>
<td class="gt_row gt_center">205 (37)</td>
<td class="gt_row gt_center">102 (36)</td>
<td class="gt_row gt_center">124 (37)</td>
<td class="gt_row gt_center">50 (39)</td>
<td class="gt_row gt_center">81 (37)</td>
<td class="gt_row gt_center">52 (34)</td></tr>
    <tr><td class="gt_row gt_left">Serotype, n (%)</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center"></td></tr>
    <tr><td class="gt_row gt_left" style="text-align: left; text-indent: 10px;">DENV-1</td>
<td class="gt_row gt_center">228 (41)</td>
<td class="gt_row gt_center">121 (43)</td>
<td class="gt_row gt_center">161 (48)</td>
<td class="gt_row gt_center">61 (48)</td>
<td class="gt_row gt_center">67 (31)</td>
<td class="gt_row gt_center">60 (39)</td></tr>
    <tr><td class="gt_row gt_left" style="text-align: left; text-indent: 10px;">DENV-2</td>
<td class="gt_row gt_center">74 (13)</td>
<td class="gt_row gt_center">47 (17)</td>
<td class="gt_row gt_center">22 (7)</td>
<td class="gt_row gt_center">16 (13)</td>
<td class="gt_row gt_center">52 (24)</td>
<td class="gt_row gt_center">31 (20)</td></tr>
    <tr><td class="gt_row gt_left" style="text-align: left; text-indent: 10px;">DENV-3</td>
<td class="gt_row gt_center">59 (11)</td>
<td class="gt_row gt_center">29 (10)</td>
<td class="gt_row gt_center">43 (13)</td>
<td class="gt_row gt_center">18 (14)</td>
<td class="gt_row gt_center">16 (7)</td>
<td class="gt_row gt_center">11 (7)</td></tr>
    <tr><td class="gt_row gt_left" style="text-align: left; text-indent: 10px;">DENV-4</td>
<td class="gt_row gt_center">161 (29)</td>
<td class="gt_row gt_center">70 (25)</td>
<td class="gt_row gt_center">91 (27)</td>
<td class="gt_row gt_center">26 (20)</td>
<td class="gt_row gt_center">70 (32)</td>
<td class="gt_row gt_center">44 (29)</td></tr>
    <tr><td class="gt_row gt_left" style="text-align: left; text-indent: 10px;">Unknown</td>
<td class="gt_row gt_center">34 (6)</td>
<td class="gt_row gt_center">14 (5)</td>
<td class="gt_row gt_center">20 (6)</td>
<td class="gt_row gt_center">6 (5)</td>
<td class="gt_row gt_center">14 (6)</td>
<td class="gt_row gt_center">8 (5)</td></tr>
    <tr><td class="gt_row gt_left">Immune status, n (%)</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center"></td></tr>
    <tr><td class="gt_row gt_left" style="text-align: left; text-indent: 10px;">Probable primary</td>
<td class="gt_row gt_center">124 (22)</td>
<td class="gt_row gt_center">41 (15)</td>
<td class="gt_row gt_center">86 (26)</td>
<td class="gt_row gt_center">15 (12)</td>
<td class="gt_row gt_center">38 (17)</td>
<td class="gt_row gt_center">26 (17)</td></tr>
    <tr><td class="gt_row gt_left" style="text-align: left; text-indent: 10px;">Probable secondary</td>
<td class="gt_row gt_center">355 (64)</td>
<td class="gt_row gt_center">218 (78)</td>
<td class="gt_row gt_center">202 (60)</td>
<td class="gt_row gt_center">100 (79)</td>
<td class="gt_row gt_center">153 (70)</td>
<td class="gt_row gt_center">118 (77)</td></tr>
    <tr><td class="gt_row gt_left" style="text-align: left; text-indent: 10px;">Inconclusive</td>
<td class="gt_row gt_center">77 (14)</td>
<td class="gt_row gt_center">22 (8)</td>
<td class="gt_row gt_center">49 (15)</td>
<td class="gt_row gt_center">12 (9)</td>
<td class="gt_row gt_center">28 (13)</td>
<td class="gt_row gt_center">10 (6)</td></tr>
    <tr><td class="gt_row gt_left">Obesity, n (%)</td>
<td class="gt_row gt_center">78 (14)</td>
<td class="gt_row gt_center">29 (10)</td>
<td class="gt_row gt_center">62 (18)</td>
<td class="gt_row gt_center">19 (15)</td>
<td class="gt_row gt_center">16 (7)</td>
<td class="gt_row gt_center">10 (6)</td></tr>
    <tr><td class="gt_row gt_left">Diabetes, n (%)</td>
<td class="gt_row gt_center">4 (1)</td>
<td class="gt_row gt_center">1 (0)</td>
<td class="gt_row gt_center">0 (0)</td>
<td class="gt_row gt_center">0 (0)</td>
<td class="gt_row gt_center">4 (2)</td>
<td class="gt_row gt_center">1 (1)</td></tr>
    <tr><td class="gt_row gt_left">WHO 2009 classification, n (%)</td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center"></td>
<td class="gt_row gt_center"></td></tr>
    <tr><td class="gt_row gt_left" style="text-align: left; text-indent: 10px;">Mild dengue</td>
<td class="gt_row gt_center">266 (48)</td>
<td class="gt_row gt_center">49 (17)</td>
<td class="gt_row gt_center">168 (50)</td>
<td class="gt_row gt_center">17 (13)</td>
<td class="gt_row gt_center">98 (45)</td>
<td class="gt_row gt_center">32 (21)</td></tr>
    <tr><td class="gt_row gt_left" style="text-align: left; text-indent: 10px;">Dengue with warning signs</td>
<td class="gt_row gt_center">288 (52)</td>
<td class="gt_row gt_center">186 (66)</td>
<td class="gt_row gt_center">169 (50)</td>
<td class="gt_row gt_center">81 (64)</td>
<td class="gt_row gt_center">119 (54)</td>
<td class="gt_row gt_center">105 (68)</td></tr>
    <tr><td class="gt_row gt_left" style="text-align: left; text-indent: 10px;">Severe dengue</td>
<td class="gt_row gt_center">0 (0)</td>
<td class="gt_row gt_center">43 (15)</td>
<td class="gt_row gt_center">0 (0)</td>
<td class="gt_row gt_center">27 (21)</td>
<td class="gt_row gt_center">0 (0)</td>
<td class="gt_row gt_center">16 (10)</td></tr>
    <tr><td class="gt_row gt_left" style="text-align: left; text-indent: 10px;">Unknown</td>
<td class="gt_row gt_center">2 (0)</td>
<td class="gt_row gt_center">3 (1)</td>
<td class="gt_row gt_center">0 (0)</td>
<td class="gt_row gt_center">2 (2)</td>
<td class="gt_row gt_center">2 (1)</td>
<td class="gt_row gt_center">1 (1)</td></tr>
    <tr><td class="gt_row gt_left">Hospitalization, n (%)</td>
<td class="gt_row gt_center">175 (31)</td>
<td class="gt_row gt_center">161 (57)</td>
<td class="gt_row gt_center">127 (38)</td>
<td class="gt_row gt_center">83 (65)</td>
<td class="gt_row gt_center">48 (22)</td>
<td class="gt_row gt_center">78 (51)</td></tr>
  </tbody>
  
  
</table>
</div>

########################################################################################### 

#### Figure 2. Biomarker levels by groups.

``` r
tick <- c(0,1,4,10,40,100,200,400,1000,4000,10000,20000,40000,70000) # for y-axis tick labels

p2 <- dat_plot %>%
  ggplot(., aes(Daygr, Result^(1/4), fill=group2, color=group2)) +
  geom_boxplot(alpha=.5, outlier.size=.9, lwd=.4, fatten=1) +
  geom_boxplot(alpha=0, outlier.color=NA, color="black", lwd=.4, fatten=1) +
  facet_wrap(~ Biomarker, scales="free", ncol=5) +
  scale_y_continuous(breaks=tick^(1/4), labels=tick) +
  xlab("Illness day (day 1 [n=140]; day 2 [n=390]; day 3 [n=307]; day 10-20 [n=625]; day >20 [n=43])") +
  theme(axis.title.y=element_blank(), legend.position="top", legend.title=element_blank(),
        axis.text.y=element_text(size=rel(.8)))
p2
```

![](Elife-ERA-codes_files/figure-gfm/fig%202-1.png)<!-- -->

``` r
#ggsave(filename="Fig 2.pdf", dpi=300, plot=p2, width=11, height=8)
```

########################################################################################### 

#### Table 2. Results from models for the primary endpoint (severe or moderate dengue).

``` r
# Set datadist for 'lrm' function [rms]
dd <- datadist(dat1); options(datadist="dd")

# Use my function 'get_est1' to get ORs and CIs from models
tmp1 <- data.frame(
  sort1 = rep(c(1:10), 2),
  sort2 = c(rep(2,10), rep(3,10)),
  bio = rep(c("VCAM", "SDC", "Ang", "IL8", "IP10", "IL1RA", "CD163", "TREM", "Fer", "CRP"), 2),
  ref1 = c(ref0-1, ref0),
  ref2 = c(ref0, ref0+1)
) %>%
  arrange(sort1, sort2) %>%
  group_by(sort1, sort2) %>%
  do(cbind(.,
           # Children - single models
           or1c = get_est1(out="sev.or.inte", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="OR", age=10, dat=dat1),
           lo1c = get_est1(out="sev.or.inte", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="loCI", age=10, dat=dat1),
           up1c = get_est1(out="sev.or.inte", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="upCI", age=10, dat=dat1),
           # Children - global model
           or2c = get_est2(out="sev.or.inte", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="OR", age=10, dat=dat1),
           lo2c = get_est2(out="sev.or.inte", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="loCI", age=10, dat=dat1),
           up2c = get_est2(out="sev.or.inte", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="upCI", age=10, dat=dat1),
           # Adults - single models
           or1a = get_est1(out="sev.or.inte", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="OR", age=25, dat=dat1),
           lo1a = get_est1(out="sev.or.inte", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="loCI", age=25, dat=dat1),
           up1a = get_est1(out="sev.or.inte", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="upCI", age=25, dat=dat1),
           # Adults - global model
           or2a = get_est2(out="sev.or.inte", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="OR", age=25, dat=dat1),
           lo2a = get_est2(out="sev.or.inte", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="loCI", age=25, dat=dat1),
           up2a = get_est2(out="sev.or.inte", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="upCI", age=25, dat=dat1))) %>%
  ungroup()
for (i in 6:17) {tmp1[[i]] <- sprintf("%.2f", round(tmp1[[i]],2))}

# Use my function 'get_est1' to get p-values from models
tmp2 <- data.frame(
  sort1 = c(1:10),
  sort2 = rep(1,10),
  bio = c("VCAM", "SDC", "Ang", "IL8", "IP10", "IL1RA", "CD163", "TREM", "Fer", "CRP")
) %>%
  group_by(sort1) %>%
  do(cbind(.,
           p.s1 = get_est1(out="sev.or.inte", bio=.$bio, est="p", dat=dat1), # P overall from single models
           p.s1_int = get_est1(out="sev.or.inte", bio=.$bio, est="p int", dat=dat1), # P interaction from single models
           p.g1 = get_est2(out="sev.or.inte", bio=.$bio, est="p", dat=dat1), # P overall from global models
           p.g1_int = get_est2(out="sev.or.inte", bio=.$bio, est="p int", dat=dat1))) %>% # P interaction from global models
  ungroup()
for (i in 4:7) {tmp2[[i]] <- ifelse(tmp2[[i]]<0.001, "<0.001", sprintf("%.3f", round(tmp2[[i]],3)))}  

# Combine results into a table
res1 <- bind_rows(tmp1, tmp2) %>%
  arrange(sort1, sort2) %>%
  mutate(bio = ifelse(!is.na(ref1), paste(" -", round(2^ref2,0), "vs", round(2^ref1,0), sep=" "), as.character(bio)),
         or.sc1 = ifelse(is.na(lo1c), NA, paste(or1c, " (", lo1c, "-", up1c, ")", sep="")), # s: single model; c: children
         or.gc1 = ifelse(is.na(lo2c), NA, paste(or2c, " (", lo2c, "-", up2c, ")", sep="")), # g: global model
         or.sa1 = ifelse(is.na(lo1a), NA, paste(or1a, " (", lo1a, "-", up1a, ")", sep="")), # a: adults
         or.ga1 = ifelse(is.na(lo2a), NA, paste(or2a, " (", lo2a, "-", up2a, ")", sep=""))) %>%
  select(bio, or.sc1, or.sa1, p.s1, p.s1_int, or.gc1, or.ga1, p.g1, p.g1_int)

names(res1) <- c("", "OR (children - single)", "OR (adults - single)", "P overall (single)", "P interaction (single)",
                 "OR (children - global)", "OR (adults - global)", "P overall (global)", "P interaction (global)")

# Report the results
knitr::kable(res1)
```

|                  | OR (children - single) | OR (adults - single) | P overall (single) | P interaction (single) | OR (children - global) | OR (adults - global) | P overall (global) | P interaction (global) |
|:-----------------|:-----------------------|:---------------------|:-------------------|:-----------------------|:-----------------------|:---------------------|:-------------------|:-----------------------|
| VCAM             |                        |                      | &lt;0.001          | 0.715                  |                        |                      | 0.441              | 0.213                  |
| \- 1636 vs 818   | 1.20 (1.04-1.38)       | 1.35 (1.15-1.58)     |                    |                        | 0.90 (0.73-1.10)       | 1.22 (0.96-1.57)     |                    |                        |
| \- 3272 vs 1636  | 1.25 (1.02-1.53)       | 1.48 (1.19-1.85)     |                    |                        | 0.87 (0.66-1.15)       | 1.30 (0.93-1.80)     |                    |                        |
| SDC              |                        |                      | &lt;0.001          | 0.088                  |                        |                      | 0.002              | 0.588                  |
| \- 2519 vs 1260  | 2.67 (1.31-5.43)       | 3.33 (1.32-8.42)     |                    |                        | 2.03 (0.77-5.34)       | 5.11 (1.56-16.78)    |                    |                        |
| \- 5039 vs 2519  | 1.71 (1.18-2.47)       | 3.71 (2.09-6.58)     |                    |                        | 1.76 (0.98-3.14)       | 2.52 (1.17-5.42)     |                    |                        |
| Ang              |                        |                      | &lt;0.001          | 0.524                  |                        |                      | 0.039              | 0.068                  |
| \- 1204 vs 602   | 1.64 (1.39-1.94)       | 1.51 (1.26-1.82)     |                    |                        | 1.67 (1.23-2.25)       | 1.01 (0.74-1.38)     |                    |                        |
| \- 2409 vs 1204  | 2.21 (1.58-3.10)       | 2.00 (1.40-2.85)     |                    |                        | 1.95 (1.25-3.05)       | 1.01 (0.65-1.57)     |                    |                        |
| IL8              |                        |                      | &lt;0.001          | &lt;0.001              |                        |                      | &lt;0.001          | &lt;0.001              |
| \- 14 vs 7       | 1.42 (1.05-1.91)       | 2.18 (1.47-3.24)     |                    |                        | 0.91 (0.63-1.34)       | 1.69 (1.05-2.71)     |                    |                        |
| \- 28 vs 14      | 0.99 (0.78-1.25)       | 2.33 (1.63-3.33)     |                    |                        | 0.53 (0.36-0.77)       | 2.05 (1.34-3.13)     |                    |                        |
| IP10             |                        |                      | &lt;0.001          | 0.984                  |                        |                      | 0.206              | 0.630                  |
| \- 3093 vs 1546  | 1.46 (1.26-1.68)       | 1.45 (1.21-1.73)     |                    |                        | 0.94 (0.73-1.19)       | 0.80 (0.57-1.12)     |                    |                        |
| \- 6186 vs 3093  | 1.68 (1.35-2.09)       | 1.69 (1.29-2.22)     |                    |                        | 1.08 (0.77-1.51)       | 0.82 (0.52-1.29)     |                    |                        |
| IL1RA            |                        |                      | &lt;0.001          | 0.082                  |                        |                      | &lt;0.001          | 0.032                  |
| \- 6434 vs 3217  | 1.69 (1.42-2.03)       | 1.48 (1.21-1.81)     |                    |                        | 2.07 (1.52-2.84)       | 1.45 (0.98-2.15)     |                    |                        |
| \- 12868 vs 6434 | 1.82 (1.46-2.27)       | 1.70 (1.29-2.24)     |                    |                        | 2.16 (1.53-3.05)       | 1.47 (0.94-2.30)     |                    |                        |
| CD163            |                        |                      | &lt;0.001          | 0.551                  |                        |                      | 0.217              | 0.341                  |
| \- 295 vs 147    | 1.57 (1.14-2.15)       | 1.49 (1.13-1.98)     |                    |                        | 1.40 (0.89-2.22)       | 1.27 (0.84-1.91)     |                    |                        |
| \- 589 vs 295    | 1.46 (1.10-1.93)       | 1.61 (1.09-2.37)     |                    |                        | 1.21 (0.87-1.69)       | 1.39 (0.89-2.18)     |                    |                        |
| TREM             |                        |                      | 0.059              | 0.997                  |                        |                      | 0.555              | 0.393                  |
| \- 85 vs 42      | 1.87 (1.23-2.84)       | 1.79 (1.10-2.93)     |                    |                        | 1.13 (0.70-1.81)       | 1.21 (0.65-2.26)     |                    |                        |
| \- 169 vs 85     | 1.12 (0.91-1.38)       | 1.12 (0.82-1.53)     |                    |                        | 0.89 (0.65-1.21)       | 0.61 (0.38-0.99)     |                    |                        |
| Fer              |                        |                      | 0.042              | 0.054                  |                        |                      | 0.008              | 0.002                  |
| \- 243 vs 122    | 1.18 (1.01-1.38)       | 1.06 (0.89-1.27)     |                    |                        | 1.30 (1.04-1.64)       | 0.78 (0.61-0.99)     |                    |                        |
| \- 487 vs 243    | 1.26 (1.00-1.58)       | 0.90 (0.66-1.23)     |                    |                        | 1.22 (0.89-1.67)       | 0.66 (0.44-1.00)     |                    |                        |
| CRP              |                        |                      | &lt;0.001          | 0.031                  |                        |                      | 0.184              | 0.138                  |
| \- 28 vs 14      | 1.26 (1.12-1.41)       | 1.25 (1.03-1.52)     |                    |                        | 1.08 (0.93-1.25)       | 1.10 (0.85-1.44)     |                    |                        |
| \- 56 vs 28      | 1.13 (0.95-1.34)       | 1.38 (1.11-1.71)     |                    |                        | 0.93 (0.75-1.15)       | 1.36 (1.02-1.81)     |                    |                        |

########################################################################################### 

#### Figure 3. Results from models for the primary endpoint (severe or moderate dengue).

``` r
# Set datadist for 'lrm' function [rms]
dd <- datadist(dat1); options(datadist="dd")

# Get results from models for children
dd$limits["Adjust to","Age"] <- 10

dat_m1 <- get_pred1(out="sev.or.inte", bio="VCAM", age=10, dat=dat1) %>% rename(value=VCAM) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_m2 <- get_pred1(out="sev.or.inte", bio="SDC", age=10, dat=dat1) %>% rename(value=SDC) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_m3 <- get_pred1(out="sev.or.inte", bio="Ang", age=10, dat=dat1) %>% rename(value=Ang) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_m4 <- get_pred1(out="sev.or.inte", bio="IL8", age=10, dat=dat1) %>% rename(value=IL8) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_m5 <- get_pred1(out="sev.or.inte", bio="IP10", age=10, dat=dat1) %>% rename(value=IP10) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_m6 <- get_pred1(out="sev.or.inte", bio="IL1RA", age=10, dat=dat1) %>% rename(value=IL1RA) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_m7 <- get_pred1(out="sev.or.inte", bio="CD163", age=10, dat=dat1) %>% rename(value=CD163) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_m8 <- get_pred1(out="sev.or.inte", bio="TREM", age=10, dat=dat1) %>% rename(value=TREM) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_m9 <- get_pred1(out="sev.or.inte", bio="Fer", age=10, dat=dat1) %>% rename(value=Fer) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_m10 <- get_pred1(out="sev.or.inte", bio="CRP", age=10, dat=dat1) %>% rename(value=CRP) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))

dat_mg1 <- get_pred2(out="sev.or.inte", bio="VCAM", age=10, dat=dat1) %>% rename(value=VCAM) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_mg2 <- get_pred2(out="sev.or.inte", bio="SDC", age=10, dat=dat1) %>% rename(value=SDC) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_mg3 <- get_pred2(out="sev.or.inte", bio="Ang", age=10, dat=dat1) %>% rename(value=Ang) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_mg4 <- get_pred2(out="sev.or.inte", bio="IL8", age=10, dat=dat1) %>% rename(value=IL8) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_mg5 <- get_pred2(out="sev.or.inte", bio="IP10", age=10, dat=dat1) %>% rename(value=IP10) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_mg6 <- get_pred2(out="sev.or.inte", bio="IL1RA", age=10, dat=dat1) %>% rename(value=IL1RA) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_mg7 <- get_pred2(out="sev.or.inte", bio="CD163", age=10, dat=dat1) %>% rename(value=CD163) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_mg8 <- get_pred2(out="sev.or.inte", bio="TREM", age=10, dat=dat1) %>% rename(value=TREM) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_mg9 <- get_pred2(out="sev.or.inte", bio="Fer", age=10, dat=dat1) %>% rename(value=Fer) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_mg10 <- get_pred2(out="sev.or.inte", bio="CRP", age=10, dat=dat1) %>% rename(value=CRP) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))

dat_p1 <- rbind(dat_m1,dat_mg1, dat_m2,dat_mg2, dat_m3,dat_mg3, dat_m4,dat_mg4, dat_m5,dat_mg5, dat_m6,dat_mg6, dat_m7,dat_mg7, dat_m8,dat_mg8, dat_m9,dat_mg9, dat_m10,dat_mg10) %>%
  mutate(age = "10 years")

# Get results from models for adults
dd$limits["Adjust to","Age"] <- 25

dat_m1 <- get_pred1(out="sev.or.inte", bio="VCAM", age=25, dat=dat1) %>% rename(value=VCAM) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_m2 <- get_pred1(out="sev.or.inte", bio="SDC", age=25, dat=dat1) %>% rename(value=SDC) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_m3 <- get_pred1(out="sev.or.inte", bio="Ang", age=25, dat=dat1) %>% rename(value=Ang) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_m4 <- get_pred1(out="sev.or.inte", bio="IL8", age=25, dat=dat1) %>% rename(value=IL8) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_m5 <- get_pred1(out="sev.or.inte", bio="IP10", age=25, dat=dat1) %>% rename(value=IP10) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_m6 <- get_pred1(out="sev.or.inte", bio="IL1RA", age=25, dat=dat1) %>% rename(value=IL1RA) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_m7 <- get_pred1(out="sev.or.inte", bio="CD163", age=25, dat=dat1) %>% rename(value=CD163) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_m8 <- get_pred1(out="sev.or.inte", bio="TREM", age=25, dat=dat1) %>% rename(value=TREM) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_m9 <- get_pred1(out="sev.or.inte", bio="Fer", age=25, dat=dat1) %>% rename(value=Fer) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_m10 <- get_pred1(out="sev.or.inte", bio="CRP", age=25, dat=dat1) %>% rename(value=CRP) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))

dat_mg1 <- get_pred2(out="sev.or.inte", bio="VCAM", age=25, dat=dat1) %>% rename(value=VCAM) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_mg2 <- get_pred2(out="sev.or.inte", bio="SDC", age=25, dat=dat1) %>% rename(value=SDC) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_mg3 <- get_pred2(out="sev.or.inte", bio="Ang", age=25, dat=dat1) %>% rename(value=Ang) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_mg4 <- get_pred2(out="sev.or.inte", bio="IL8", age=25, dat=dat1) %>% rename(value=IL8) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_mg5 <- get_pred2(out="sev.or.inte", bio="IP10", age=25, dat=dat1) %>% rename(value=IP10) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_mg6 <- get_pred2(out="sev.or.inte", bio="IL1RA", age=25, dat=dat1) %>% rename(value=IL1RA) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_mg7 <- get_pred2(out="sev.or.inte", bio="CD163", age=25, dat=dat1) %>% rename(value=CD163) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_mg8 <- get_pred2(out="sev.or.inte", bio="TREM", age=25, dat=dat1) %>% rename(value=TREM) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_mg9 <- get_pred2(out="sev.or.inte", bio="Fer", age=25, dat=dat1) %>% rename(value=Fer) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_mg10 <- get_pred2(out="sev.or.inte", bio="CRP", age=25, dat=dat1) %>% rename(value=CRP) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))

dat_p2 <- rbind(dat_m1,dat_mg1, dat_m2,dat_mg2, dat_m3,dat_mg3, dat_m4,dat_mg4, dat_m5,dat_mg5, dat_m6,dat_mg6, dat_m7,dat_mg7, dat_m8,dat_mg8, dat_m9,dat_mg9, dat_m10,dat_mg10) %>%
  mutate(age = "25 years")

# Merge results for children and adults
vis1 <- rbind(dat_p1, dat_p2) %>% mutate(outcome = "outcome 1")

# Merge data for plots
tmp0 <- dat1 %>% arrange(age15) %>% select(sev.or.inte)
tmp <- data.frame(sev.or.inte = rep(tmp0$sev.or.inte, 20))
  
tmp_p1 <- vis1 %>%
  arrange(biomarker, model) %>%
  mutate(value1 = 2^value,
         age = factor(age, levels=c("10 years", "25 years")),
         gr = ifelse(model=="Single model" & age=="10 years", 1,
                     ifelse(model=="Single model" & age=="25 years", 2,
                            ifelse(model=="Global model" & age=="10 years", 3, 4))),
         gr = factor(gr, levels=c(1:4), labels=c("Single children", "Single adults", "Global children", "Global adults")),
         model = factor(model, levels=c("Single model", "Global model"))) %>%
  bind_cols(., tmp)

vline <- tmp_p1 %>%
  filter(model=="Single model") %>%
  filter(yhat==0) %>%
  select(biomarker, value, value1) %>%
  rename(vline=value, vline1=value1)

dat_p <- left_join(tmp_p1, vline, by="biomarker")

# Alternative data to limit to 5th - 95th of the values
dat_alt <- dat_p %>%
  group_by(biomarker, model) %>%
  mutate(up = quantile(value, .95),
         lo = quantile(value, .05)) %>%
  ungroup() %>%
  mutate(is.outlier = value<lo | value>up | value<0,
         value = ifelse(is.outlier, NA, value),
         value1 = ifelse(is.outlier, NA, value1),
         yhat = ifelse(is.outlier, NA, yhat),
         lower = ifelse(is.outlier, NA, lower),
         upper = ifelse(is.outlier, NA, upper))

dat_1_5 <- dat_alt %>% 
  filter(biomarker %in% c("VCAM", "SDC", "Ang", "IL8", "IP10")) %>%
  mutate(biomarker = factor(biomarker, levels=c("VCAM", "SDC", "Ang", "IL8", "IP10")))

dat_6_10 <- dat_alt %>% 
  filter(biomarker %in% c("IL1RA", "CD163", "TREM", "Fer", "CRP")) %>%
  mutate(biomarker = factor(biomarker, levels=c("IL1RA", "CD163", "TREM", "Fer", "CRP")))

# Modify facets' scales
#require(facetscales)
xVCAM <- c(1,4,15,60,250,1000,4000)
xSDC <- c(1400,2000,2800,4000,5600)
xAng <- c(50,100,250,500,1000,2000)
xIL8 <- c(5,7,10,14,20,28,40)
xIP10 <- c(25,100,400,1600,6400)
xIL1RA <- c(1000,2000,4000,8000,16000)
xCD163 <- c(75,150,300,600)
xTREM <- c(35,50,70,100,140,200)
xFer <- c(50,100,200,400,800)
xCRP <- c(2.5,5,10,20,40,80)

scales_x <- list(
  `VCAM` = scale_x_continuous(breaks = log2(xVCAM), labels = xVCAM),
  `SDC` = scale_x_continuous(breaks = log2(xSDC), labels = xSDC),
  `Ang` = scale_x_continuous(breaks = log2(xAng), labels = xAng), 
  `IL8` = scale_x_continuous(breaks = log2(xIL8), labels = xIL8), 
  `IP10` = scale_x_continuous(breaks = log2(xIP10), labels = xIP10), 
  `IL1RA` = scale_x_continuous(breaks = log2(xIL1RA), labels = xIL1RA), 
  `CD163` = scale_x_continuous(breaks = log2(xCD163), labels = xCD163), 
  `TREM` = scale_x_continuous(breaks = log2(xTREM), labels = xTREM), 
  `Fer` = scale_x_continuous(breaks = log2(xFer), labels = xFer), 
  `CRP` = scale_x_continuous(breaks = log2(xCRP), labels = xCRP)
)

ybreak1 <- c(.125, .25, .5, 1, 2, 4, 8)
ybreak2 <- c(.125, .25, .5, 1, 2, 4)

scales_y1 <- list(
  `Single model` = scale_y_continuous(limits = c(log(.124), log(8.01)), breaks = log(ybreak1), labels = ybreak1),
  `Global model` = scale_y_continuous(limits = c(log(.124), log(8.01)), breaks = log(ybreak1), labels = ybreak1)
)

scales_y2 <- list(
  `Single model` = scale_y_continuous(limits = c(log(.124), log(4.01)), breaks = log(ybreak2), labels = ybreak2),
  `Global model` = scale_y_continuous(limits = c(log(.124), log(4.01)), breaks = log(ybreak2), labels = ybreak2)
)

# Set facets' names
models <- c(`Single model` = "Single model",
            `Global model` = "Global model")

biomarkers <- c(`VCAM` = "VCAM-1 (ng/ml)",
                `SDC`= "SDC-1 (pg/ml)",
                `Ang` = "Ang-2 (pg/ml)",
                `IL8` = "IL-8 (pg/ml)",
                `IP10` = "IP-10 (pg/ml)",
                `IL1RA` = "IL-1RA (pg/ml)",
                `CD163` = "sCD163 (ng/ml)",
                `TREM` = "sTREM-1 (pg/ml)",
                `Fer` = "Ferritin (ng/ml)",
                `CRP` = "CRP (mg/l)")

# Plot for the first 5 biomarkers
p.1.5 <- ggplot(dat_1_5, aes(x=value)) +
  geom_ribbon(aes(ymin=lower, ymax=upper, fill=age), alpha=.15) +
  geom_vline(aes(xintercept=vline), color="black", linetype="dashed", alpha=.4) +
  geom_hline(yintercept=0, color="black", linetype="dashed", alpha=.4) +
  geom_line(aes(y=yhat, color=age), size=.5, alpha=.7) +
  geom_rug(data=filter(dat_1_5, model=="Single model" & sev.or.inte==0 & age=="10 years"), sides="b", alpha=.2, color="red") +
  geom_rug(data=filter(dat_1_5, model=="Single model" & sev.or.inte==1 & age=="10 years"), sides="t", alpha=.2, color="red") +
  geom_rug(data=filter(dat_1_5, model=="Global model" & sev.or.inte==0 & age=="25 years"), sides="b", alpha=.2, color="blue") +
  geom_rug(data=filter(dat_1_5, model=="Global model" & sev.or.inte==1 & age=="25 years"), sides="t", alpha=.2, color="blue") +
  scale_color_manual(values=c("red","blue"), labels=c("Children","Adults")) +
  scale_fill_manual(values=c("red","blue"), labels=c("Children","Adults")) +
  ylab("Odds ratio") +
  theme(legend.position="top", legend.title=element_blank(), axis.title.x=element_blank(), axis.text.x=element_text(angle=45)) +
  facet_grid_sc(cols=vars(biomarker), rows=vars(model), 
                scales=list(x=scales_x, y=scales_y1),
                labeller = as_labeller(c(models, biomarkers)))

# Plot for the last 5 biomarkers
p.6.10 <- ggplot(dat_6_10, aes(x=value)) +
  geom_ribbon(aes(ymin=lower, ymax=upper, fill=age), alpha=.15) +
  geom_vline(aes(xintercept=vline), color="black", linetype="dashed", alpha=.4) +
  geom_hline(yintercept=0, color="black", linetype="dashed", alpha=.4) +
  geom_line(aes(y=yhat, color=age), size=.5, alpha=.7) +
  geom_rug(data=filter(dat_6_10, model=="Single model" & sev.or.inte==0 & age=="10 years"), sides="b", alpha=.2, color="red") +
  geom_rug(data=filter(dat_6_10, model=="Single model" & sev.or.inte==1 & age=="10 years"), sides="t", alpha=.2, color="red") +
  geom_rug(data=filter(dat_6_10, model=="Global model" & sev.or.inte==0 & age=="25 years"), sides="b", alpha=.2, color="blue") +
  geom_rug(data=filter(dat_6_10, model=="Global model" & sev.or.inte==1 & age=="25 years"), sides="t", alpha=.2, color="blue") +
  scale_color_manual(values=c("red","blue"), labels=c("Children","Adults")) +
  scale_fill_manual(values=c("red","blue"), labels=c("Children","Adults")) +
  ylab("Odds ratio") +
  theme(legend.position="none", axis.title.x=element_blank(), axis.text.x=element_text(angle=45)) +
  facet_grid_sc(cols=vars(biomarker), rows=vars(model), 
                scales=list(x=scales_x, y=scales_y2),
                labeller = as_labeller(c(models, biomarkers)))

# Merge plots
p3 <- gridExtra::grid.arrange(p.1.5, p.6.10, nrow=2, heights=c(1,.86))
```

![](Elife-ERA-codes_files/figure-gfm/fig%203-1.png)<!-- -->

``` r
#ggsave(filename="Fig 3.pdf", p3, dpi=300, width=11, height=8)
```

########################################################################################### 

#### Table 3. Best combinations of biomarkers associated with severe or moderate dengue for children.

``` r
# EPV --------------------------------------------------------------
pred <- c("VCAM", "SDC", "Ang", "IL8", "IP10", "IL1RA", "CD163", "TREM", "Fer", "CRP")

# Estimate full model ----------------------------------------------
full_mod <- glm(sev.or.inte ~ VCAM + SDC + Ang + IL8 + IP10 + IL1RA + CD163 + TREM + Fer + CRP, 
                family=binomial, data=dat1c, x=T, y=T)

# Selected model ---------------------------------------------------
sel_var <- matrix(0, ncol=length(pred)+1, nrow=5, dimnames=list(NULL, c(pred, "aic")))

for (i in 1:5) {
  if (i==1) {
    bs <- dredge(full_mod, rank="AIC")
  } else {
    bs <- dredge(full_mod, rank="AIC", m.lim=c(i,i))
  }
  bs_var <- attr(get.models(bs, 1)[[1]]$terms, "term.labels")

  for (j in 1:(ncol(sel_var)-1)) {sel_var[i,j] <- ifelse(names(sel_var[i,j]) %in% bs_var, 1, 0)}
  formula <- paste("sev.or.inte~", paste(names(sel_var[i,][sel_var[i,]==1]), collapse = "+"))
  sel_mod <- glm(formula, data = dat1c, family = binomial, x = T, y = T)
  sel_var[i, ncol(sel_var)] <- AIC(sel_mod)
}

# Report results --------------------------------------------------
out1 <- as.data.frame(sel_var) %>% mutate(AIC = round(aic,1)) %>% select(-aic)
for (i in 1:(ncol(out1)-1)) {out1[,i] <- ifelse(out1[,i]==0, NA, "+")}

out2 <- as.data.frame(t(out1))
colnames(out2) <- c("Best of all combinations", "Best combination of 2 variables", "Best combination of 3 variables",
                    "Best combination of 4 variables", "Best combination of 5 variables")
rownames(out2) <- c("- VCAM-1", "- SDC-1", "- Ang-2", "- IL-8", "- IP-10", "- IL-1RA", "- sCD163", "- sTREM-1",
                    "- Ferritin", "- CRP", "AIC of the selected model")

knitr::kable(out2)
```

|                           | Best of all combinations | Best combination of 2 variables | Best combination of 3 variables | Best combination of 4 variables | Best combination of 5 variables |
|:--------------------------|:-------------------------|:--------------------------------|:--------------------------------|:--------------------------------|:--------------------------------|
| \- VCAM-1                 |                          |                                 |                                 |                                 |                                 |
| \- SDC-1                  | \+                       |                                 |                                 |                                 |                                 |
| \- Ang-2                  | \+                       |                                 | \+                              | \+                              | \+                              |
| \- IL-8                   | \+                       |                                 |                                 |                                 | \+                              |
| \- IP-10                  | \+                       |                                 |                                 | \+                              | \+                              |
| \- IL-1RA                 | \+                       | \+                              | \+                              | \+                              | \+                              |
| \- sCD163                 |                          |                                 |                                 |                                 |                                 |
| \- sTREM-1                |                          |                                 |                                 |                                 |                                 |
| \- Ferritin               | \+                       | \+                              | \+                              | \+                              | \+                              |
| \- CRP                    |                          |                                 |                                 |                                 |                                 |
| AIC of the selected model | 465.9                    | 484.7                           | 480.0                           | 473.7                           | 467.6                           |

``` r
# For bootstrap results please look at Appendix 7-tables 1, 2 (for the best of all combinations)
# and the bootstrap codes for the best combinations of 2, 3, 4, and 5 variables
```

########################################################################################### 

#### Table 4. Best combinations of biomarkers associated with severe or moderate dengue for adults.

``` r
# EPV --------------------------------------------------------------
pred <- c("VCAM", "SDC", "Ang", "IL8", "ns1(IP10)", "IL1RA", "CD163", "TREM", "Fer", "CRP")

# Estimate full model ----------------------------------------------
full_mod <- glm(sev.or.inte ~ VCAM + SDC + Ang + IL8 + ns1(IP10) + IL1RA + CD163 + TREM + Fer + CRP, 
                family=binomial, data=dat1a, x=T, y=T)

# Selected model ---------------------------------------------------
sel_var <- matrix(0, ncol=length(pred)+1, nrow=5, dimnames=list(NULL, c(pred, "aic")))

for (i in 1:5) {
  if (i==1) {
    bs <- dredge(full_mod, rank="AIC")
  } else {
    bs <- dredge(full_mod, rank="AIC", m.lim=c(i,i))
  }
  bs_var <- attr(get.models(bs, 1)[[1]]$terms, "term.labels")

  for (j in 1:(ncol(sel_var)-1)) {sel_var[i,j] <- ifelse(names(sel_var[i,j]) %in% bs_var, 1, 0)}
  formula <- paste("sev.or.inte~", paste(names(sel_var[i,][sel_var[i,]==1]), collapse = "+"))
  sel_mod <- glm(formula, data = dat1a, family = binomial, x = T, y = T)
  sel_var[i, ncol(sel_var)] <- AIC(sel_mod)
}

# Report results --------------------------------------------------
out1 <- as.data.frame(sel_var) %>% mutate(AIC = round(aic,1)) %>% select(-aic)
for (i in 1:(ncol(out1)-1)) {out1[,i] <- ifelse(out1[,i]==0, NA, "+")}

out2 <- as.data.frame(t(out1))
colnames(out2) <- c("Best of all combinations", "Best combination of 2 variables", "Best combination of 3 variables",
                    "Best combination of 4 variables", "Best combination of 5 variables")
rownames(out2) <- c("- VCAM-1", "- SDC-1", "- Ang-2", "- IL-8", "- IP-10", "- IL-1RA", "- sCD163", "- sTREM-1",
                    "- Ferritin", "- CRP", "AIC of the selected model")

knitr::kable(out2)
```

|                           | Best of all combinations | Best combination of 2 variables | Best combination of 3 variables | Best combination of 4 variables | Best combination of 5 variables |
|:--------------------------|:-------------------------|:--------------------------------|:--------------------------------|:--------------------------------|:--------------------------------|
| \- VCAM-1                 |                          |                                 |                                 |                                 |                                 |
| \- SDC-1                  | \+                       | \+                              | \+                              | \+                              | \+                              |
| \- Ang-2                  |                          |                                 |                                 |                                 |                                 |
| \- IL-8                   | \+                       | \+                              | \+                              | \+                              | \+                              |
| \- IP-10                  | \+                       |                                 |                                 |                                 |                                 |
| \- IL-1RA                 | \+                       |                                 |                                 | \+                              | \+                              |
| \- sCD163                 | \+                       |                                 |                                 |                                 |                                 |
| \- sTREM-1                | \+                       |                                 |                                 |                                 | \+                              |
| \- Ferritin               | \+                       |                                 | \+                              | \+                              | \+                              |
| \- CRP                    |                          |                                 |                                 |                                 |                                 |
| AIC of the selected model | 430.5                    | 441.1                           | 434.2                           | 431.6                           | 430.7                           |

``` r
# For bootstrap results please look at Appendix 7-tables 3, 4 (for the best of all combinations)
# and the bootstrap codes for the best combinations of 2, 3, 4, and 5 variables
```

########################################################################################### 

#### Appendix 3. Statistical analysis. Analysis to find the best combination of biomarkers to predict the primary endpoint (aim \#2) - Step \#1 - AICs of models

``` r
# For children
screen_child <- data.frame(bio = c("VCAM", "SDC", "Ang", "IL8", "IP10", "IL1RA", "CD163", "TREM", "Fer", "CRP")) %>%
  group_by(bio) %>%
  do(cbind(.,
           model1_child = screen(mod=1, bio=.$bio, dat=dat1c),
           model2_child = screen(mod=2, bio=.$bio, dat=dat1c),
           model3_child = screen(mod=3, bio=.$bio, dat=dat1c),
           model4_child = screen(mod=4, bio=.$bio, dat=dat1c))) %>%
  ungroup() %>%
  mutate(bio = factor(bio, levels=c("VCAM", "SDC", "Ang", "IL8", "IP10", "IL1RA", "CD163", "TREM", "Fer", "CRP"))) %>%
  arrange(bio)

# For adults
screen_adult <- data.frame(bio = c("VCAM", "SDC", "Ang", "IL8", "IP10", "IL1RA", "CD163", "TREM", "Fer", "CRP")) %>%
  group_by(bio) %>%
  do(cbind(.,
           model1_adult = screen(mod=1, bio=.$bio, dat=dat1a),
           model2_adult = screen(mod=2, bio=.$bio, dat=dat1a),
           model3_adult = screen(mod=3, bio=.$bio, dat=dat1a),
           model4_adult = screen(mod=4, bio=.$bio, dat=dat1a))) %>%
  ungroup() %>%
  mutate(bio = factor(bio, levels=c("VCAM", "SDC", "Ang", "IL8", "IP10", "IL1RA", "CD163", "TREM", "Fer", "CRP"))) %>%
  arrange(bio)

# Combine results
out <- full_join(screen_child, screen_adult, by="bio")
knitr::kable(out)
```

| bio   | model1\_child | model2\_child | model3\_child | model4\_child | model1\_adult | model2\_adult | model3\_adult | model4\_adult |
|:------|--------------:|--------------:|--------------:|--------------:|--------------:|--------------:|--------------:|--------------:|
| VCAM  |         530.8 |         531.3 |         530.6 |         532.6 |         499.8 |         496.7 |         501.3 |         496.4 |
| SDC   |         537.9 |         538.0 |               |               |         459.5 |         461.0 |               |               |
| Ang   |         511.0 |         509.6 |         512.5 |         509.7 |         493.5 |         490.8 |         493.3 |         492.8 |
| IL8   |         548.5 |         549.0 |         548.1 |         546.8 |         457.7 |         457.8 |         456.5 |         457.8 |
| IP10  |         521.2 |         517.2 |         523.0 |         518.3 |         500.7 |         492.6 |         502.5 |         494.5 |
| IL1RA |         492.9 |         494.9 |               |               |         493.5 |         494.2 |               |               |
| CD163 |         531.9 |         533.2 |         531.9 |         533.2 |         497.6 |         499.3 |         499.5 |         501.1 |
| TREM  |         545.7 |         544.2 |         545.6 |         545.4 |         507.8 |         509.5 |         507.0 |         505.0 |
| Fer   |         542.6 |         544.4 |               |               |         509.3 |         505.2 |               |               |
| CRP   |         536.7 |         536.4 |               |               |         505.1 |         506.8 |               |               |

########################################################################################### 

#### Appendix 4—table 1. Summary of clinical phenotype of the primary endpoint.

``` r
dat %>%
  filter(sev.or.inte == 1) %>%
  mutate(age15 = factor(age15, levels = c("No", "Yes"), labels = c("Children", "Adults")),
         mod.only = sev.or.inte - sev.only) %>% # to create 'moderate dengue'
  select(age15, sev.only, sev.leak, shock, res.dis, sev.neu, sev.bleed, sev.other, sev.liver, mod.only, mod.leak, mod.liver, mod.bleed, mod.other, mod.neu) %>%
  tbl_summary(by = age15,
              statistic = list(all_categorical() ~ "{n} ({p})"),
              label = list(sev.only ~ "Severe dengue, n (%)", sev.leak ~ "- Severe plasma leakage", 
                           shock ~ " + Dengue shock syndrome", res.dis ~ " + Respiratory distress", 
                           sev.neu ~ "- Severe neurologic involvement", sev.bleed ~ "- Severe bleeding", 
                           sev.other ~ "- Severe other major organ failure", sev.liver ~ "- Severe hepatic involvement", 
                           mod.only ~ "Moderate dengue, n (%)", mod.leak ~ "- Moderate plasma leakage", 
                           mod.liver ~ "- Moderate hepatic involvement", mod.bleed ~ "- Moderate bleeding", 
                           mod.other ~ "- Moderate other major organ involvement", mod.neu ~ "- Moderate neurologic involvement")) %>%
  add_overall() %>%
  modify_header(label = "", stat_0 = "**All patients (N={N})**", stat_by = "**{level} (N={n})**") %>%
  modify_footnote(update = everything() ~ NA)
```

<div id="akjangducp" style="overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>html {
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Helvetica Neue', 'Fira Sans', 'Droid Sans', Arial, sans-serif;
}

#akjangducp .gt_table {
  display: table;
  border-collapse: collapse;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 16px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: auto;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}

#akjangducp .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#akjangducp .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#akjangducp .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 0;
  padding-bottom: 4px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#akjangducp .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#akjangducp .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#akjangducp .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}

#akjangducp .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}

#akjangducp .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#akjangducp .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#akjangducp .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}

#akjangducp .gt_group_heading {
  padding: 8px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
}

#akjangducp .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}

#akjangducp .gt_from_md > :first-child {
  margin-top: 0;
}

#akjangducp .gt_from_md > :last-child {
  margin-bottom: 0;
}

#akjangducp .gt_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}

#akjangducp .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 12px;
}

#akjangducp .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#akjangducp .gt_first_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
}

#akjangducp .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#akjangducp .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#akjangducp .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#akjangducp .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#akjangducp .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#akjangducp .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding: 4px;
}

#akjangducp .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#akjangducp .gt_sourcenote {
  font-size: 90%;
  padding: 4px;
}

#akjangducp .gt_left {
  text-align: left;
}

#akjangducp .gt_center {
  text-align: center;
}

#akjangducp .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#akjangducp .gt_font_normal {
  font-weight: normal;
}

#akjangducp .gt_font_bold {
  font-weight: bold;
}

#akjangducp .gt_font_italic {
  font-style: italic;
}

#akjangducp .gt_super {
  font-size: 65%;
}

#akjangducp .gt_footnote_marks {
  font-style: italic;
  font-weight: normal;
  font-size: 65%;
}
</style>
<table class="gt_table">
  
  <thead class="gt_col_headings">
    <tr>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="1" colspan="1"></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1"><strong>All patients (N=281)</strong></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1"><strong>Children (N=127)</strong></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1"><strong>Adults (N=154)</strong></th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr><td class="gt_row gt_left">Severe dengue, n (%)</td>
<td class="gt_row gt_center">38 (14)</td>
<td class="gt_row gt_center">29 (23)</td>
<td class="gt_row gt_center">9 (6)</td></tr>
    <tr><td class="gt_row gt_left">- Severe plasma leakage</td>
<td class="gt_row gt_center">33 (12)</td>
<td class="gt_row gt_center">24 (19)</td>
<td class="gt_row gt_center">9 (6)</td></tr>
    <tr><td class="gt_row gt_left"> + Dengue shock syndrome</td>
<td class="gt_row gt_center">25 (9)</td>
<td class="gt_row gt_center">18 (14)</td>
<td class="gt_row gt_center">7 (5)</td></tr>
    <tr><td class="gt_row gt_left"> + Respiratory distress</td>
<td class="gt_row gt_center">12 (4)</td>
<td class="gt_row gt_center">9 (7)</td>
<td class="gt_row gt_center">3 (2)</td></tr>
    <tr><td class="gt_row gt_left">- Severe neurologic involvement</td>
<td class="gt_row gt_center">3 (1)</td>
<td class="gt_row gt_center">3 (2)</td>
<td class="gt_row gt_center">0 (0)</td></tr>
    <tr><td class="gt_row gt_left">- Severe bleeding</td>
<td class="gt_row gt_center">2 (1)</td>
<td class="gt_row gt_center">2 (2)</td>
<td class="gt_row gt_center">0 (0)</td></tr>
    <tr><td class="gt_row gt_left">- Severe other major organ failure</td>
<td class="gt_row gt_center">1 (0)</td>
<td class="gt_row gt_center">0 (0)</td>
<td class="gt_row gt_center">1 (1)</td></tr>
    <tr><td class="gt_row gt_left">- Severe hepatic involvement</td>
<td class="gt_row gt_center">0 (0)</td>
<td class="gt_row gt_center">0 (0)</td>
<td class="gt_row gt_center">0 (0)</td></tr>
    <tr><td class="gt_row gt_left">Moderate dengue, n (%)</td>
<td class="gt_row gt_center">243 (86)</td>
<td class="gt_row gt_center">98 (77)</td>
<td class="gt_row gt_center">145 (94)</td></tr>
    <tr><td class="gt_row gt_left">- Moderate plasma leakage</td>
<td class="gt_row gt_center">159 (57)</td>
<td class="gt_row gt_center">73 (57)</td>
<td class="gt_row gt_center">86 (56)</td></tr>
    <tr><td class="gt_row gt_left">- Moderate hepatic involvement</td>
<td class="gt_row gt_center">102 (36)</td>
<td class="gt_row gt_center">35 (28)</td>
<td class="gt_row gt_center">67 (44)</td></tr>
    <tr><td class="gt_row gt_left">- Moderate bleeding</td>
<td class="gt_row gt_center">9 (3)</td>
<td class="gt_row gt_center">3 (2)</td>
<td class="gt_row gt_center">6 (4)</td></tr>
    <tr><td class="gt_row gt_left">- Moderate other major organ involvement</td>
<td class="gt_row gt_center">1 (0)</td>
<td class="gt_row gt_center">0 (0)</td>
<td class="gt_row gt_center">1 (1)</td></tr>
    <tr><td class="gt_row gt_left">- Moderate neurologic involvement</td>
<td class="gt_row gt_center">0 (0)</td>
<td class="gt_row gt_center">0 (0)</td>
<td class="gt_row gt_center">0 (0)</td></tr>
  </tbody>
  
  
</table>
</div>

########################################################################################### 

#### Appendix 4—table 2. Summary of biomarkers’ data (at enrollment)

``` r
# All patients
t1 <- dat0 %>%
  filter(Time == "Enrolment") %>%
  mutate(viremia = Viremia/(10^6)) %>% # to summarize Viremia as 10^6 copies/ml
  select(group2, VCAM1, SDC1, Ang2, IL8, IP10, IL1RA, CD163, TREM1, Fer, CRP, viremia) %>%
  tbl_summary(by = group2,
              statistic = list(all_continuous() ~ "{median} ({p25}, {p75})"),
              digits = list(all_continuous() ~ c(0,0), viremia ~ c(1,1)),
              label = list(VCAM1 ~ "VCAM-1 (ng/ml)", SDC1 ~ "SDC-1 (pg/ml)", Ang2 ~ "Ang-2 (pg/ml)", IL8 ~ "IL-8 (pg/ml)",
                           IP10 ~ "IP-10 (pg/ml)", IL1RA ~ "IL-1RA (pg/ml)", CD163 ~ "sCD163 (ng/ml)", TREM1 ~ "sTREM-1 (pg/ml)",
                           Fer ~ "Ferritin (ng/ml)", CRP ~ "CRP (mg/l)", viremia ~ "Viremia (10^6 copies/ml)")) %>%
  modify_header(label = "", stat_by = "**{level} (N={n})**") %>%
  modify_footnote(update = everything() ~ NA)

# Children (<15 years of age)
t2 <- dat0 %>%
  filter(Time == "Enrolment") %>%
  filter(age15 == "No") %>%
  mutate(viremia = Viremia/(10^6)) %>% # to summarize Viremia as 10^6 copies/ml
  select(group2, VCAM1, SDC1, Ang2, IL8, IP10, IL1RA, CD163, TREM1, Fer, CRP, viremia) %>%
  tbl_summary(by = group2,
              statistic = list(all_continuous() ~ "{median} ({p25}, {p75})"),
              digits = list(all_continuous() ~ c(0,0), viremia ~ c(1,1)),
              label = list(VCAM1 ~ "VCAM-1 (ng/ml)", SDC1 ~ "SDC-1 (pg/ml)", Ang2 ~ "Ang-2 (pg/ml)", IL8 ~ "IL-8 (pg/ml)",
                           IP10 ~ "IP-10 (pg/ml)", IL1RA ~ "IL-1RA (pg/ml)", CD163 ~ "sCD163 (ng/ml)", TREM1 ~ "sTREM-1 (pg/ml)",
                           Fer ~ "Ferritin (ng/ml)", CRP ~ "CRP (mg/l)", viremia ~ "Viremia (10^6 copies/ml)")) %>%
  modify_header(label = "", stat_by = "**{level} (N={n})**") %>%
  modify_footnote(update = everything() ~ NA)

# Adults (>=15 years of age)
t3 <- dat0 %>%
  filter(Time == "Enrolment") %>%
  filter(age15 == "Yes") %>%
  mutate(viremia = Viremia/(10^6)) %>% # to summarize Viremia as 10^6 copies/ml
  select(group2, VCAM1, SDC1, Ang2, IL8, IP10, IL1RA, CD163, TREM1, Fer, CRP, viremia) %>%
  tbl_summary(by = group2,
              statistic = list(all_continuous() ~ "{median} ({p25}, {p75})"),
              digits = list(all_continuous() ~ c(0,0), viremia ~ c(1,1)),
              label = list(VCAM1 ~ "VCAM-1 (ng/ml)", SDC1 ~ "SDC-1 (pg/ml)", Ang2 ~ "Ang-2 (pg/ml)", IL8 ~ "IL-8 (pg/ml)",
                           IP10 ~ "IP-10 (pg/ml)", IL1RA ~ "IL-1RA (pg/ml)", CD163 ~ "sCD163 (ng/ml)", TREM1 ~ "sTREM-1 (pg/ml)",
                           Fer ~ "Ferritin (ng/ml)", CRP ~ "CRP (mg/l)", viremia ~ "Viremia (10^6 copies/ml)")) %>%
  modify_header(label = "", stat_by = "**{level} (N={n})**") %>%
  modify_footnote(update = everything() ~ NA)

# Merge tables
tbl_merge(tbls = list(t1, t2, t3),
          tab_spanner = c("**All patients**", "**Children**", "**Adults**"))
```

<div id="rchxvpddux" style="overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>html {
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Helvetica Neue', 'Fira Sans', 'Droid Sans', Arial, sans-serif;
}

#rchxvpddux .gt_table {
  display: table;
  border-collapse: collapse;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 16px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: auto;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}

#rchxvpddux .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#rchxvpddux .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#rchxvpddux .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 0;
  padding-bottom: 4px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#rchxvpddux .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#rchxvpddux .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#rchxvpddux .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}

#rchxvpddux .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}

#rchxvpddux .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#rchxvpddux .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#rchxvpddux .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}

#rchxvpddux .gt_group_heading {
  padding: 8px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
}

#rchxvpddux .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}

#rchxvpddux .gt_from_md > :first-child {
  margin-top: 0;
}

#rchxvpddux .gt_from_md > :last-child {
  margin-bottom: 0;
}

#rchxvpddux .gt_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}

#rchxvpddux .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 12px;
}

#rchxvpddux .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#rchxvpddux .gt_first_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
}

#rchxvpddux .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#rchxvpddux .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#rchxvpddux .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#rchxvpddux .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#rchxvpddux .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#rchxvpddux .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding: 4px;
}

#rchxvpddux .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#rchxvpddux .gt_sourcenote {
  font-size: 90%;
  padding: 4px;
}

#rchxvpddux .gt_left {
  text-align: left;
}

#rchxvpddux .gt_center {
  text-align: center;
}

#rchxvpddux .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#rchxvpddux .gt_font_normal {
  font-weight: normal;
}

#rchxvpddux .gt_font_bold {
  font-weight: bold;
}

#rchxvpddux .gt_font_italic {
  font-style: italic;
}

#rchxvpddux .gt_super {
  font-size: 65%;
}

#rchxvpddux .gt_footnote_marks {
  font-style: italic;
  font-weight: normal;
  font-size: 65%;
}
</style>
<table class="gt_table">
  
  <thead class="gt_col_headings">
    <tr>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="2" colspan="1"></th>
      <th class="gt_center gt_columns_top_border gt_column_spanner_outer" rowspan="1" colspan="2">
        <span class="gt_column_spanner"><strong>All patients</strong></span>
      </th>
      <th class="gt_center gt_columns_top_border gt_column_spanner_outer" rowspan="1" colspan="2">
        <span class="gt_column_spanner"><strong>Children</strong></span>
      </th>
      <th class="gt_center gt_columns_top_border gt_column_spanner_outer" rowspan="1" colspan="2">
        <span class="gt_column_spanner"><strong>Adults</strong></span>
      </th>
    </tr>
    <tr>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1"><strong>Uncomplicated dengue (N=556)</strong></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1"><strong>Severe/moderate dengue (N=281)</strong></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1"><strong>Uncomplicated dengue (N=337)</strong></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1"><strong>Severe/moderate dengue (N=127)</strong></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1"><strong>Uncomplicated dengue (N=219)</strong></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1"><strong>Severe/moderate dengue (N=154)</strong></th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr><td class="gt_row gt_left">VCAM-1 (ng/ml)</td>
<td class="gt_row gt_center">1,404 (540, 2,548)</td>
<td class="gt_row gt_center">2,027 (1,122, 3,577)</td>
<td class="gt_row gt_center">1,442 (447, 2,546)</td>
<td class="gt_row gt_center">2,020 (1,232, 3,384)</td>
<td class="gt_row gt_center">1,356 (568, 2,560)</td>
<td class="gt_row gt_center">2,092 (1,060, 4,202)</td></tr>
    <tr><td class="gt_row gt_left">SDC-1 (pg/ml)</td>
<td class="gt_row gt_center">2,334 (1,864, 3,131)</td>
<td class="gt_row gt_center">2,997 (2,230, 4,201)</td>
<td class="gt_row gt_center">2,369 (1,861, 3,423)</td>
<td class="gt_row gt_center">2,846 (2,164, 4,173)</td>
<td class="gt_row gt_center">2,260 (1,879, 2,898)</td>
<td class="gt_row gt_center">3,122 (2,278, 4,211)</td></tr>
    <tr><td class="gt_row gt_left">Ang-2 (pg/ml)</td>
<td class="gt_row gt_center">1,064 (550, 1,584)</td>
<td class="gt_row gt_center">1,521 (899, 2,318)</td>
<td class="gt_row gt_center">1,102 (584, 1,563)</td>
<td class="gt_row gt_center">1,547 (967, 2,318)</td>
<td class="gt_row gt_center">944 (516, 1,585)</td>
<td class="gt_row gt_center">1,516 (885, 2,321)</td></tr>
    <tr><td class="gt_row gt_left">IL-8 (pg/ml)</td>
<td class="gt_row gt_center">12 (8, 22)</td>
<td class="gt_row gt_center">17 (11, 28)</td>
<td class="gt_row gt_center">15 (9, 26)</td>
<td class="gt_row gt_center">16 (10, 27)</td>
<td class="gt_row gt_center">10 (7, 15)</td>
<td class="gt_row gt_center">19 (12, 29)</td></tr>
    <tr><td class="gt_row gt_left">IP-10 (pg/ml)</td>
<td class="gt_row gt_center">2,502 (732, 4,509)</td>
<td class="gt_row gt_center">4,092 (2,436, 6,441)</td>
<td class="gt_row gt_center">2,245 (458, 4,531)</td>
<td class="gt_row gt_center">3,942 (2,046, 6,287)</td>
<td class="gt_row gt_center">2,793 (1,370, 4,495)</td>
<td class="gt_row gt_center">4,242 (2,524, 6,469)</td></tr>
    <tr><td class="gt_row gt_left">IL-1RA (pg/ml)</td>
<td class="gt_row gt_center">5,237 (2,603, 9,082)</td>
<td class="gt_row gt_center">9,105 (5,933, 14,977)</td>
<td class="gt_row gt_center">4,491 (2,318, 8,977)</td>
<td class="gt_row gt_center">9,688 (6,109, 16,786)</td>
<td class="gt_row gt_center">5,721 (3,479, 9,703)</td>
<td class="gt_row gt_center">8,993 (5,953, 12,935)</td></tr>
    <tr><td class="gt_row gt_left">sCD163 (ng/ml)</td>
<td class="gt_row gt_center">278 (185, 447)</td>
<td class="gt_row gt_center">322 (228, 503)</td>
<td class="gt_row gt_center">326 (212, 481)</td>
<td class="gt_row gt_center">386 (256, 603)</td>
<td class="gt_row gt_center">226 (157, 374)</td>
<td class="gt_row gt_center">291 (207, 410)</td></tr>
    <tr><td class="gt_row gt_left">sTREM-1 (pg/ml)</td>
<td class="gt_row gt_center">81 (59, 114)</td>
<td class="gt_row gt_center">96 (69, 132)</td>
<td class="gt_row gt_center">80 (58, 115)</td>
<td class="gt_row gt_center">93 (67, 128)</td>
<td class="gt_row gt_center">84 (60, 114)</td>
<td class="gt_row gt_center">98 (73, 134)</td></tr>
    <tr><td class="gt_row gt_left">Ferritin (ng/ml)</td>
<td class="gt_row gt_center">233 (116, 406)</td>
<td class="gt_row gt_center">261 (133, 433)</td>
<td class="gt_row gt_center">177 (99, 324)</td>
<td class="gt_row gt_center">224 (110, 402)</td>
<td class="gt_row gt_center">303 (161, 510)</td>
<td class="gt_row gt_center">278 (160, 448)</td></tr>
    <tr><td class="gt_row gt_left">CRP (mg/l)</td>
<td class="gt_row gt_center">25 (10, 54)</td>
<td class="gt_row gt_center">34 (17, 72)</td>
<td class="gt_row gt_center">18 (7, 41)</td>
<td class="gt_row gt_center">24 (13, 58)</td>
<td class="gt_row gt_center">38 (17, 65)</td>
<td class="gt_row gt_center">45 (25, 80)</td></tr>
    <tr><td class="gt_row gt_left">Viremia (10^6 copies/ml)</td>
<td class="gt_row gt_center">15.8 (0.7, 148.5)</td>
<td class="gt_row gt_center">79.2 (5.3, 582.7)</td>
<td class="gt_row gt_center">21.8 (1.8, 167.0)</td>
<td class="gt_row gt_center">105.4 (8.4, 646.0)</td>
<td class="gt_row gt_center">9.8 (0.3, 115.5)</td>
<td class="gt_row gt_center">56.2 (3.6, 496.0)</td></tr>
  </tbody>
  
  
</table>
</div>

#### Appendix 4—table 2. Summary of biomarkers’ data (at follow-up)

``` r
# All patients
t1 <- dat0 %>%
  filter(Time == "Follow up") %>%
  select(group2, VCAM1, SDC1, Ang2, IL8, IP10, IL1RA, CD163, TREM1, Fer, CRP) %>%
  tbl_summary(by = group2,
              statistic = list(all_continuous() ~ "{median} ({p25}, {p75})"),
              digits = list(all_continuous() ~ c(0,0), c(IL8, CRP) ~ c(1,1)),
              label = list(VCAM1 ~ "VCAM-1 (ng/ml)", SDC1 ~ "SDC-1 (pg/ml)", Ang2 ~ "Ang-2 (pg/ml)", IL8 ~ "IL-8 (pg/ml)",
                           IP10 ~ "IP-10 (pg/ml)", IL1RA ~ "IL-1RA (pg/ml)", CD163 ~ "sCD163 (ng/ml)", TREM1 ~ "sTREM-1 (pg/ml)",
                           Fer ~ "Ferritin (ng/ml)", CRP ~ "CRP (mg/l)")) %>%
  modify_header(label = "", stat_by = "**{level} (N={n})**") %>%
  modify_footnote(update = everything() ~ NA)

# Children (<15 years of age)
t2 <- dat0 %>%
  filter(Time == "Follow up") %>%
  filter(age15 == "No") %>%
  select(group2, VCAM1, SDC1, Ang2, IL8, IP10, IL1RA, CD163, TREM1, Fer, CRP) %>%
  tbl_summary(by = group2,
              statistic = list(all_continuous() ~ "{median} ({p25}, {p75})"),
              digits = list(all_continuous() ~ c(0,0), c(IL8, CRP) ~ c(1,1)),
              label = list(VCAM1 ~ "VCAM-1 (ng/ml)", SDC1 ~ "SDC-1 (pg/ml)", Ang2 ~ "Ang-2 (pg/ml)", IL8 ~ "IL-8 (pg/ml)",
                           IP10 ~ "IP-10 (pg/ml)", IL1RA ~ "IL-1RA (pg/ml)", CD163 ~ "sCD163 (ng/ml)", TREM1 ~ "sTREM-1 (pg/ml)",
                           Fer ~ "Ferritin (ng/ml)", CRP ~ "CRP (mg/l)")) %>%
  modify_header(label = "", stat_by = "**{level} (N={n})**") %>%
  modify_footnote(update = everything() ~ NA)

# Adults (>=15 years of age)
t3 <- dat0 %>%
  filter(Time == "Follow up") %>%
  filter(age15 == "Yes") %>%
  select(group2, VCAM1, SDC1, Ang2, IL8, IP10, IL1RA, CD163, TREM1, Fer, CRP) %>%
  tbl_summary(by = group2,
              statistic = list(all_continuous() ~ "{median} ({p25}, {p75})"),
              digits = list(all_continuous() ~ c(0,0), c(IL8, CRP) ~ c(1,1)),
              label = list(VCAM1 ~ "VCAM-1 (ng/ml)", SDC1 ~ "SDC-1 (pg/ml)", Ang2 ~ "Ang-2 (pg/ml)", IL8 ~ "IL-8 (pg/ml)",
                           IP10 ~ "IP-10 (pg/ml)", IL1RA ~ "IL-1RA (pg/ml)", CD163 ~ "sCD163 (ng/ml)", TREM1 ~ "sTREM-1 (pg/ml)",
                           Fer ~ "Ferritin (ng/ml)", CRP ~ "CRP (mg/l)")) %>%
  modify_header(label = "", stat_by = "**{level} (N={n})**") %>%
  modify_footnote(update = everything() ~ NA)

# Merge tables
tbl_merge(tbls = list(t1, t2, t3),
          tab_spanner = c("**All patients**", "**Children**", "**Adults**"))
```

<div id="jmmhwmlbth" style="overflow-x:auto;overflow-y:auto;width:auto;height:auto;">
<style>html {
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Helvetica Neue', 'Fira Sans', 'Droid Sans', Arial, sans-serif;
}

#jmmhwmlbth .gt_table {
  display: table;
  border-collapse: collapse;
  margin-left: auto;
  margin-right: auto;
  color: #333333;
  font-size: 16px;
  font-weight: normal;
  font-style: normal;
  background-color: #FFFFFF;
  width: auto;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #A8A8A8;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #A8A8A8;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
}

#jmmhwmlbth .gt_heading {
  background-color: #FFFFFF;
  text-align: center;
  border-bottom-color: #FFFFFF;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#jmmhwmlbth .gt_title {
  color: #333333;
  font-size: 125%;
  font-weight: initial;
  padding-top: 4px;
  padding-bottom: 4px;
  border-bottom-color: #FFFFFF;
  border-bottom-width: 0;
}

#jmmhwmlbth .gt_subtitle {
  color: #333333;
  font-size: 85%;
  font-weight: initial;
  padding-top: 0;
  padding-bottom: 4px;
  border-top-color: #FFFFFF;
  border-top-width: 0;
}

#jmmhwmlbth .gt_bottom_border {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#jmmhwmlbth .gt_col_headings {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
}

#jmmhwmlbth .gt_col_heading {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  padding-left: 5px;
  padding-right: 5px;
  overflow-x: hidden;
}

#jmmhwmlbth .gt_column_spanner_outer {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: normal;
  text-transform: inherit;
  padding-top: 0;
  padding-bottom: 0;
  padding-left: 4px;
  padding-right: 4px;
}

#jmmhwmlbth .gt_column_spanner_outer:first-child {
  padding-left: 0;
}

#jmmhwmlbth .gt_column_spanner_outer:last-child {
  padding-right: 0;
}

#jmmhwmlbth .gt_column_spanner {
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: bottom;
  padding-top: 5px;
  padding-bottom: 6px;
  overflow-x: hidden;
  display: inline-block;
  width: 100%;
}

#jmmhwmlbth .gt_group_heading {
  padding: 8px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
}

#jmmhwmlbth .gt_empty_group_heading {
  padding: 0.5px;
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  vertical-align: middle;
}

#jmmhwmlbth .gt_from_md > :first-child {
  margin-top: 0;
}

#jmmhwmlbth .gt_from_md > :last-child {
  margin-bottom: 0;
}

#jmmhwmlbth .gt_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  margin: 10px;
  border-top-style: solid;
  border-top-width: 1px;
  border-top-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 1px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 1px;
  border-right-color: #D3D3D3;
  vertical-align: middle;
  overflow-x: hidden;
}

#jmmhwmlbth .gt_stub {
  color: #333333;
  background-color: #FFFFFF;
  font-size: 100%;
  font-weight: initial;
  text-transform: inherit;
  border-right-style: solid;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
  padding-left: 12px;
}

#jmmhwmlbth .gt_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#jmmhwmlbth .gt_first_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
}

#jmmhwmlbth .gt_grand_summary_row {
  color: #333333;
  background-color: #FFFFFF;
  text-transform: inherit;
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
}

#jmmhwmlbth .gt_first_grand_summary_row {
  padding-top: 8px;
  padding-bottom: 8px;
  padding-left: 5px;
  padding-right: 5px;
  border-top-style: double;
  border-top-width: 6px;
  border-top-color: #D3D3D3;
}

#jmmhwmlbth .gt_striped {
  background-color: rgba(128, 128, 128, 0.05);
}

#jmmhwmlbth .gt_table_body {
  border-top-style: solid;
  border-top-width: 2px;
  border-top-color: #D3D3D3;
  border-bottom-style: solid;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
}

#jmmhwmlbth .gt_footnotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#jmmhwmlbth .gt_footnote {
  margin: 0px;
  font-size: 90%;
  padding: 4px;
}

#jmmhwmlbth .gt_sourcenotes {
  color: #333333;
  background-color: #FFFFFF;
  border-bottom-style: none;
  border-bottom-width: 2px;
  border-bottom-color: #D3D3D3;
  border-left-style: none;
  border-left-width: 2px;
  border-left-color: #D3D3D3;
  border-right-style: none;
  border-right-width: 2px;
  border-right-color: #D3D3D3;
}

#jmmhwmlbth .gt_sourcenote {
  font-size: 90%;
  padding: 4px;
}

#jmmhwmlbth .gt_left {
  text-align: left;
}

#jmmhwmlbth .gt_center {
  text-align: center;
}

#jmmhwmlbth .gt_right {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

#jmmhwmlbth .gt_font_normal {
  font-weight: normal;
}

#jmmhwmlbth .gt_font_bold {
  font-weight: bold;
}

#jmmhwmlbth .gt_font_italic {
  font-style: italic;
}

#jmmhwmlbth .gt_super {
  font-size: 65%;
}

#jmmhwmlbth .gt_footnote_marks {
  font-style: italic;
  font-weight: normal;
  font-size: 65%;
}
</style>
<table class="gt_table">
  
  <thead class="gt_col_headings">
    <tr>
      <th class="gt_col_heading gt_columns_bottom_border gt_left" rowspan="2" colspan="1"></th>
      <th class="gt_center gt_columns_top_border gt_column_spanner_outer" rowspan="1" colspan="2">
        <span class="gt_column_spanner"><strong>All patients</strong></span>
      </th>
      <th class="gt_center gt_columns_top_border gt_column_spanner_outer" rowspan="1" colspan="2">
        <span class="gt_column_spanner"><strong>Children</strong></span>
      </th>
      <th class="gt_center gt_columns_top_border gt_column_spanner_outer" rowspan="1" colspan="2">
        <span class="gt_column_spanner"><strong>Adults</strong></span>
      </th>
    </tr>
    <tr>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1"><strong>Uncomplicated dengue (N=437)</strong></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1"><strong>Severe/moderate dengue (N=231)</strong></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1"><strong>Uncomplicated dengue (N=292)</strong></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1"><strong>Severe/moderate dengue (N=112)</strong></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1"><strong>Uncomplicated dengue (N=145)</strong></th>
      <th class="gt_col_heading gt_columns_bottom_border gt_center" rowspan="1" colspan="1"><strong>Severe/moderate dengue (N=119)</strong></th>
    </tr>
  </thead>
  <tbody class="gt_table_body">
    <tr><td class="gt_row gt_left">VCAM-1 (ng/ml)</td>
<td class="gt_row gt_center">402 (102, 730)</td>
<td class="gt_row gt_center">686 (344, 961)</td>
<td class="gt_row gt_center">579 (182, 858)</td>
<td class="gt_row gt_center">782 (402, 1,078)</td>
<td class="gt_row gt_center">173 (26, 388)</td>
<td class="gt_row gt_center">622 (343, 835)</td></tr>
    <tr><td class="gt_row gt_left">SDC-1 (pg/ml)</td>
<td class="gt_row gt_center">2,769 (2,298, 3,514)</td>
<td class="gt_row gt_center">3,417 (2,815, 5,495)</td>
<td class="gt_row gt_center">2,957 (2,319, 4,115)</td>
<td class="gt_row gt_center">3,122 (2,748, 5,507)</td>
<td class="gt_row gt_center">2,666 (2,196, 3,058)</td>
<td class="gt_row gt_center">3,745 (2,971, 5,495)</td></tr>
    <tr><td class="gt_row gt_left">Ang-2 (pg/ml)</td>
<td class="gt_row gt_center">953 (478, 1,479)</td>
<td class="gt_row gt_center">1,155 (675, 1,567)</td>
<td class="gt_row gt_center">1,163 (738, 1,646)</td>
<td class="gt_row gt_center">1,352 (710, 1,856)</td>
<td class="gt_row gt_center">565 (302, 923)</td>
<td class="gt_row gt_center">1,044 (626, 1,345)</td></tr>
    <tr><td class="gt_row gt_left">IL-8 (pg/ml)</td>
<td class="gt_row gt_center">4.9 (2.3, 12.4)</td>
<td class="gt_row gt_center">5.7 (2.7, 10.4)</td>
<td class="gt_row gt_center">6.8 (3.0, 15.1)</td>
<td class="gt_row gt_center">5.9 (2.4, 10.5)</td>
<td class="gt_row gt_center">2.7 (1.6, 4.8)</td>
<td class="gt_row gt_center">5.5 (3.1, 10.4)</td></tr>
    <tr><td class="gt_row gt_left">IP-10 (pg/ml)</td>
<td class="gt_row gt_center">57 (24, 91)</td>
<td class="gt_row gt_center">76 (47, 133)</td>
<td class="gt_row gt_center">67 (33, 98)</td>
<td class="gt_row gt_center">86 (38, 143)</td>
<td class="gt_row gt_center">39 (22, 70)</td>
<td class="gt_row gt_center">75 (48, 108)</td></tr>
    <tr><td class="gt_row gt_left">IL-1RA (pg/ml)</td>
<td class="gt_row gt_center">412 (279, 635)</td>
<td class="gt_row gt_center">455 (328, 626)</td>
<td class="gt_row gt_center">441 (323, 687)</td>
<td class="gt_row gt_center">501 (352, 664)</td>
<td class="gt_row gt_center">336 (210, 480)</td>
<td class="gt_row gt_center">407 (308, 615)</td></tr>
    <tr><td class="gt_row gt_left">sCD163 (ng/ml)</td>
<td class="gt_row gt_center">337 (216, 553)</td>
<td class="gt_row gt_center">412 (257, 661)</td>
<td class="gt_row gt_center">340 (226, 562)</td>
<td class="gt_row gt_center">456 (279, 680)</td>
<td class="gt_row gt_center">328 (199, 523)</td>
<td class="gt_row gt_center">386 (241, 589)</td></tr>
    <tr><td class="gt_row gt_left">sTREM-1 (pg/ml)</td>
<td class="gt_row gt_center">99 (73, 132)</td>
<td class="gt_row gt_center">90 (68, 116)</td>
<td class="gt_row gt_center">98 (72, 132)</td>
<td class="gt_row gt_center">91 (67, 115)</td>
<td class="gt_row gt_center">101 (73, 134)</td>
<td class="gt_row gt_center">90 (70, 116)</td></tr>
    <tr><td class="gt_row gt_left">Ferritin (ng/ml)</td>
<td class="gt_row gt_center">202 (120, 309)</td>
<td class="gt_row gt_center">273 (181, 382)</td>
<td class="gt_row gt_center">177 (112, 263)</td>
<td class="gt_row gt_center">209 (154, 311)</td>
<td class="gt_row gt_center">267 (160, 404)</td>
<td class="gt_row gt_center">322 (247, 436)</td></tr>
    <tr><td class="gt_row gt_left">CRP (mg/l)</td>
<td class="gt_row gt_center">0.7 (0.3, 1.8)</td>
<td class="gt_row gt_center">0.8 (0.4, 2.0)</td>
<td class="gt_row gt_center">0.6 (0.3, 1.3)</td>
<td class="gt_row gt_center">0.6 (0.3, 1.1)</td>
<td class="gt_row gt_center">1.1 (0.5, 2.7)</td>
<td class="gt_row gt_center">1.1 (0.5, 3.4)</td></tr>
    <tr><td class="gt_row gt_left" style="text-align: left; text-indent: 10px;">Unknown</td>
<td class="gt_row gt_center">1</td>
<td class="gt_row gt_center">3</td>
<td class="gt_row gt_center">0</td>
<td class="gt_row gt_center">1</td>
<td class="gt_row gt_center">1</td>
<td class="gt_row gt_center">2</td></tr>
  </tbody>
  
  
</table>
</div>

########################################################################################### 

#### Appendix 4-figure 1. Biomarker levels by individual.

``` r
tick <- c(0,1,10,40,100,200,400,1000,4000,10000,20000,40000,70000) # for y-axis tick labels

p41 <- ggplot(dat_plot, aes(Day, Result^(1/4), color=group2)) +
  geom_point(alpha = 0.5, size = .5) +
  geom_line(aes(group = Code), alpha = 0.07) +
  facet_wrap(~ Biomarker, scales="free", ncol=5) +
  scale_y_continuous(breaks=tick^(1/4), labels=tick) +
  scale_x_continuous(breaks = c(1,3,10,15,20,25,30), 
                     limit = c(1,31),
                     name = "Illness day") +
  theme_bw() +
  theme(axis.title.y=element_blank(), legend.position="top", legend.title=element_blank(),
        axis.text.y=element_text(size=rel(.8)))
p41
```

![](Elife-ERA-codes_files/figure-gfm/fig%20A4.1-1.png)<!-- -->

``` r
#ggsave(filename="A4_fig1.pdf", plot=p41, width=8.5, height=6.2)
```

########################################################################################### 

#### Appendix 4—figure 2. Pairwise correlation of biomarker levels at enrollment and age.

``` r
collabel <- c("Age", "VCAM-1", "SDC-1", "Ang-2", "IL-8", "IP-10", "IL-1RA", "sCD163", "sTREM-1", "Ferritin", "CRP", "Viremia")

# Change biomarker's values to log-2 and viremia to log
dat_fig42 <- dat %>%
  mutate(log2_VCAM1 = as.numeric(log2(VCAM1)),
         log2_SDC1 = as.numeric(log2(SDC1)),
         log2_Ang2 = as.numeric(log2(Ang2)),
         log2_IL8 = as.numeric(log2(IL8)),
         log2_IP10 = as.numeric(log2(IP10)),
         log2_IL1RA = as.numeric(log2(IL1RA)),
         log2_CD163 = as.numeric(log2(CD163)),
         log2_TREM1 = as.numeric(log2(TREM1)),
         log2_Fer = as.numeric(log2(Fer)),
         log2_CRP = as.numeric(log2(CRP)),
         log_vir = as.numeric(log10(Viremia))) %>%
  select(Age, log2_VCAM1, log2_SDC1, log2_Ang2, log2_IL8, log2_IP10, log2_IL1RA, log2_CD163, log2_TREM1, log2_Fer, log2_CRP, log_vir) %>%
  as.data.frame(.)

p42 <- GGally::ggpairs(dat_fig42,
             lower = list(continuous = GGscatterPlot), # GGscatterPlot is a customized function
             upper = "blank",
             columnLabels = collabel) +
  theme(panel.grid.minor = element_blank())
p42
```

![](Elife-ERA-codes_files/figure-gfm/fig%20A4.2-1.png)<!-- -->

``` r
#ggsave(filename="A4_fig2.pdf", plot=p42, width=9.5, height=9.5)
```

########################################################################################### 

#### Appendix 5—figure 1. Results from single models for severe/moderate dengue with the interaction with serotype.

``` r
# Set datadist for 'lrm' function [rms]
dd <- datadist(dat1); options(datadist="dd")

# For DENV-1
dd$limits["Adjust to","Serotype"] <- "DENV-1"

## For children
dd$limits["Adjust to","Age"] <- 10
dat_m1 <- get_pred1s(out="sev.or.inte", bio="VCAM", age=10, serotype="DENV-1", dat=dat1) %>% rename(value=VCAM) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='DENV-1')) %>% mutate(age="Children")
dat_m2 <- get_pred1s(out="sev.or.inte", bio="SDC", age=10, serotype="DENV-1", dat=dat1) %>% rename(value=SDC) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='DENV-1')) %>% mutate(age="Children")
dat_m3 <- get_pred1s(out="sev.or.inte", bio="Ang", age=10, serotype="DENV-1", dat=dat1) %>% rename(value=Ang) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='DENV-1')) %>% mutate(age="Children")
dat_m4 <- get_pred1s(out="sev.or.inte", bio="IL8", age=10, serotype="DENV-1", dat=dat1) %>% rename(value=IL8) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='DENV-1')) %>% mutate(age="Children")
dat_m5 <- get_pred1s(out="sev.or.inte", bio="IP10", age=10, serotype="DENV-1", dat=dat1) %>% rename(value=IP10) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='DENV-1')) %>% mutate(age="Children")
dat_m6 <- get_pred1s(out="sev.or.inte", bio="IL1RA", age=10, serotype="DENV-1", dat=dat1) %>% rename(value=IL1RA) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='DENV-1')) %>% mutate(age="Children")
dat_m7 <- get_pred1s(out="sev.or.inte", bio="CD163", age=10, serotype="DENV-1", dat=dat1) %>% rename(value=CD163) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='DENV-1')) %>% mutate(age="Children")
dat_m8 <- get_pred1s(out="sev.or.inte", bio="TREM", age=10, serotype="DENV-1", dat=dat1) %>% rename(value=TREM) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='DENV-1')) %>% mutate(age="Children")
dat_m9 <- get_pred1s(out="sev.or.inte", bio="Fer", age=10, serotype="DENV-1", dat=dat1) %>% rename(value=Fer) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='DENV-1')) %>% mutate(age="Children")
dat_m10 <- get_pred1s(out="sev.or.inte", bio="CRP", age=10, serotype="DENV-1", dat=dat1) %>% rename(value=CRP) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='DENV-1')) %>% mutate(age="Children")

## For adults
dd$limits["Adjust to","Age"] <- 25
dat_mg1 <- get_pred1s(out="sev.or.inte", bio="VCAM", age=25, serotype="DENV-1", dat=dat1) %>% rename(value=VCAM) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='DENV-1')) %>% mutate(age="Adults")
dat_mg2 <- get_pred1s(out="sev.or.inte", bio="SDC", age=25, serotype="DENV-1", dat=dat1) %>% rename(value=SDC) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='DENV-1')) %>% mutate(age="Adults")
dat_mg3 <- get_pred1s(out="sev.or.inte", bio="Ang", age=25, serotype="DENV-1", dat=dat1) %>% rename(value=Ang) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='DENV-1')) %>% mutate(age="Adults")
dat_mg4 <- get_pred1s(out="sev.or.inte", bio="IL8", age=25, serotype="DENV-1", dat=dat1) %>% rename(value=IL8) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='DENV-1')) %>% mutate(age="Adults")
dat_mg5 <- get_pred1s(out="sev.or.inte", bio="IP10", age=25, serotype="DENV-1", dat=dat1) %>% rename(value=IP10) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='DENV-1')) %>% mutate(age="Adults")
dat_mg6 <- get_pred1s(out="sev.or.inte", bio="IL1RA", age=25, serotype="DENV-1", dat=dat1) %>% rename(value=IL1RA) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='DENV-1')) %>% mutate(age="Adults")
dat_mg7 <- get_pred1s(out="sev.or.inte", bio="CD163", age=25, serotype="DENV-1", dat=dat1) %>% rename(value=CD163) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='DENV-1')) %>% mutate(age="Adults")
dat_mg8 <- get_pred1s(out="sev.or.inte", bio="TREM", age=25, serotype="DENV-1", dat=dat1) %>% rename(value=TREM) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='DENV-1')) %>% mutate(age="Adults")
dat_mg9 <- get_pred1s(out="sev.or.inte", bio="Fer", age=25, serotype="DENV-1", dat=dat1) %>% rename(value=Fer) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='DENV-1')) %>% mutate(age="Adults")
dat_mg10 <- get_pred1s(out="sev.or.inte", bio="CRP", age=25, serotype="DENV-1", dat=dat1) %>% rename(value=CRP) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='DENV-1')) %>% mutate(age="Adults")

dat_p1 <- rbind(dat_m1,dat_mg1, dat_m2,dat_mg2, dat_m3,dat_mg3, dat_m4,dat_mg4, dat_m5,dat_mg5, dat_m6,dat_mg6, dat_m7,dat_mg7, dat_m8,dat_mg8, dat_m9,dat_mg9, dat_m10,dat_mg10) %>%
  mutate(serotype = "DENV-1")

# For others
dd$limits["Adjust to","Serotype"] <- "Others"

## For children
dd$limits["Adjust to","Age"] <- 10
dat_m1 <- get_pred1s(out="sev.or.inte", bio="VCAM", age=10, serotype="Others", dat=dat1) %>% rename(value=VCAM) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='Others')) %>% mutate(age="Children")
dat_m2 <- get_pred1s(out="sev.or.inte", bio="SDC", age=10, serotype="Others", dat=dat1) %>% rename(value=SDC) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='Others')) %>% mutate(age="Children")
dat_m3 <- get_pred1s(out="sev.or.inte", bio="Ang", age=10, serotype="Others", dat=dat1) %>% rename(value=Ang) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='Others')) %>% mutate(age="Children")
dat_m4 <- get_pred1s(out="sev.or.inte", bio="IL8", age=10, serotype="Others", dat=dat1) %>% rename(value=IL8) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='Others')) %>% mutate(age="Children")
dat_m5 <- get_pred1s(out="sev.or.inte", bio="IP10", age=10, serotype="Others", dat=dat1) %>% rename(value=IP10) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='Others')) %>% mutate(age="Children")
dat_m6 <- get_pred1s(out="sev.or.inte", bio="IL1RA", age=10, serotype="Others", dat=dat1) %>% rename(value=IL1RA) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='Others')) %>% mutate(age="Children")
dat_m7 <- get_pred1s(out="sev.or.inte", bio="CD163", age=10, serotype="Others", dat=dat1) %>% rename(value=CD163) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='Others')) %>% mutate(age="Children")
dat_m8 <- get_pred1s(out="sev.or.inte", bio="TREM", age=10, serotype="Others", dat=dat1) %>% rename(value=TREM) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='Others')) %>% mutate(age="Children")
dat_m9 <- get_pred1s(out="sev.or.inte", bio="Fer", age=10, serotype="Others", dat=dat1) %>% rename(value=Fer) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='Others')) %>% mutate(age="Children")
dat_m10 <- get_pred1s(out="sev.or.inte", bio="CRP", age=10, serotype="Others", dat=dat1) %>% rename(value=CRP) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='Others')) %>% mutate(age="Children")

## For adults
dd$limits["Adjust to","Age"] <- 25
dat_mg1 <- get_pred1s(out="sev.or.inte", bio="VCAM", age=25, serotype="Others", dat=dat1) %>% rename(value=VCAM) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='Others')) %>% mutate(age="Adults")
dat_mg2 <- get_pred1s(out="sev.or.inte", bio="SDC", age=25, serotype="Others", dat=dat1) %>% rename(value=SDC) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='Others')) %>% mutate(age="Adults")
dat_mg3 <- get_pred1s(out="sev.or.inte", bio="Ang", age=25, serotype="Others", dat=dat1) %>% rename(value=Ang) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='Others')) %>% mutate(age="Adults")
dat_mg4 <- get_pred1s(out="sev.or.inte", bio="IL8", age=25, serotype="Others", dat=dat1) %>% rename(value=IL8) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='Others')) %>% mutate(age="Adults")
dat_mg5 <- get_pred1s(out="sev.or.inte", bio="IP10", age=25, serotype="Others", dat=dat1) %>% rename(value=IP10) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='Others')) %>% mutate(age="Adults")
dat_mg6 <- get_pred1s(out="sev.or.inte", bio="IL1RA", age=25, serotype="Others", dat=dat1) %>% rename(value=IL1RA) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='Others')) %>% mutate(age="Adults")
dat_mg7 <- get_pred1s(out="sev.or.inte", bio="CD163", age=25, serotype="Others", dat=dat1) %>% rename(value=CD163) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='Others')) %>% mutate(age="Adults")
dat_mg8 <- get_pred1s(out="sev.or.inte", bio="TREM", age=25, serotype="Others", dat=dat1) %>% rename(value=TREM) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='Others')) %>% mutate(age="Adults")
dat_mg9 <- get_pred1s(out="sev.or.inte", bio="Fer", age=25, serotype="Others", dat=dat1) %>% rename(value=Fer) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='Others')) %>% mutate(age="Adults")
dat_mg10 <- get_pred1s(out="sev.or.inte", bio="CRP", age=25, serotype="Others", dat=dat1) %>% rename(value=CRP) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='Others')) %>% mutate(age="Adults")

dat_p2 <- rbind(dat_m1,dat_mg1, dat_m2,dat_mg2, dat_m3,dat_mg3, dat_m4,dat_mg4, dat_m5,dat_mg5, dat_m6,dat_mg6, dat_m7,dat_mg7, dat_m8,dat_mg8, dat_m9,dat_mg9, dat_m10,dat_mg10) %>%
  mutate(serotype = "Others")

# Merge data for children and adults
vis1 <- rbind(dat_p1, dat_p2) %>% mutate(outcome = "outcome 1")

# Merge data for plots
tmp0 <- dat1 %>% filter(!is.na(Serotype)) %>% arrange(Serotype) %>% select(sev.or.inte)
tmp <- data.frame(sev.or.inte = rep(tmp0$sev.or.inte, 20))
  
tmp_p1 <- vis1 %>%
  arrange(biomarker, age) %>%
  mutate(value1 = 2^value,
         serotype = factor(serotype, levels=c("DENV-1", "Others")),
         age = factor(age, levels=c("Children", "Adults"))) %>%
  bind_cols(., tmp)

vline <- tmp_p1 %>%
  filter(yhat==0) %>%
  select(biomarker, value, value1) %>%
  rename(vline=value, vline1=value1)

dat_p <- left_join(tmp_p1, vline, by="biomarker")

# Alternative data to limit to 5th - 95th of the values
dat_alt <- dat_p %>%
  group_by(biomarker, model) %>%
  mutate(up = quantile(value, .95),
         lo = quantile(value, .05)) %>%
  ungroup() %>%
  mutate(is.outlier = value<lo | value>up | value<0,
         value = ifelse(is.outlier, NA, value),
         value1 = ifelse(is.outlier, NA, value1),
         yhat = ifelse(is.outlier, NA, yhat),
         lower = ifelse(is.outlier, NA, lower),
         upper = ifelse(is.outlier, NA, upper))

dat_1_5 <- dat_alt %>% 
  filter(biomarker %in% c("VCAM", "SDC", "Ang", "IL8", "IP10")) %>%
  mutate(biomarker = factor(biomarker, levels=c("VCAM", "SDC", "Ang", "IL8", "IP10")))

dat_6_10 <- dat_alt %>% 
  filter(biomarker %in% c("IL1RA", "CD163", "TREM", "Fer", "CRP")) %>%
  mutate(biomarker = factor(biomarker, levels=c("IL1RA", "CD163", "TREM", "Fer", "CRP")))

# Modify facets' scales
#require(facetscales)

xVCAM <- c(1,4,15,60,250,1000,4000)
xSDC <- c(1400,2000,2800,4000,5600)
xAng <- c(50,100,250,500,1000,2000)
xIL8 <- c(5,7,10,14,20,28,40)
xIP10 <- c(25,100,400,1600,6400)
xIL1RA <- c(1000,2000,4000,8000,16000)
xCD163 <- c(75,150,300,600)
xTREM <- c(35,50,70,100,140,200)
xFer <- c(50,100,200,400,800)
xCRP <- c(2.5,5,10,20,40,80)

scales_x <- list(
  `VCAM` = scale_x_continuous(breaks = log2(xVCAM), labels = xVCAM),
  `SDC` = scale_x_continuous(breaks = log2(xSDC), labels = xSDC),
  `Ang` = scale_x_continuous(breaks = log2(xAng), labels = xAng), 
  `IL8` = scale_x_continuous(breaks = log2(xIL8), labels = xIL8), 
  `IP10` = scale_x_continuous(breaks = log2(xIP10), labels = xIP10), 
  `IL1RA` = scale_x_continuous(breaks = log2(xIL1RA), labels = xIL1RA), 
  `CD163` = scale_x_continuous(breaks = log2(xCD163), labels = xCD163), 
  `TREM` = scale_x_continuous(breaks = log2(xTREM), labels = xTREM), 
  `Fer` = scale_x_continuous(breaks = log2(xFer), labels = xFer), 
  `CRP` = scale_x_continuous(breaks = log2(xCRP), labels = xCRP)
)

ybreak1 <- c(.125, .25, .5, 1, 2, 4, 8)
ybreak2 <- c(.125, .25, .5, 1, 2, 4)

scales_y1 <- list(
  `Children` = scale_y_continuous(limits = c(log(.124), log(8.01)), breaks = log(ybreak1), labels = ybreak1),
  `Adults` = scale_y_continuous(limits = c(log(.124), log(8.01)), breaks = log(ybreak1), labels = ybreak1)
)

scales_y2 <- list(
  `Children` = scale_y_continuous(limits = c(log(.124), log(4.01)), breaks = log(ybreak2), labels = ybreak2),
  `Adults` = scale_y_continuous(limits = c(log(.124), log(4.01)), breaks = log(ybreak2), labels = ybreak2)
)

# Set facets' names
age <- c(`Children` = "Children",
         `Adults` = "Adults")

biomarkers <- c(`VCAM` = "VCAM-1 (ng/ml)",
                `SDC`= "SDC-1 (pg/ml)",
                `Ang` = "Ang-2 (pg/ml)",
                `IL8` = "IL-8 (pg/ml)",
                `IP10` = "IP-10 (pg/ml)",
                `IL1RA` = "IL-1RA (pg/ml)",
                `CD163` = "sCD163 (ng/ml)",
                `TREM` = "sTREM-1 (pg/ml)",
                `Fer` = "Ferritin (ng/ml)",
                `CRP` = "CRP (mg/l)")

# Plot for the first 5 biomarkers
p.1.5 <- ggplot(dat_1_5, aes(x=value)) +
  geom_ribbon(aes(ymin=lower, ymax=upper, fill=serotype), alpha=.15) +
  geom_vline(aes(xintercept=vline), color="black", linetype="dashed", alpha=.4) +
  geom_hline(yintercept=0, color="black", linetype="dashed", alpha=.4) +
  geom_line(aes(y=yhat, color=serotype), size=.5, alpha=.7) +
  geom_rug(data=filter(dat_1_5, age=="Children" & sev.or.inte==0 & serotype=="DENV-1"), sides="b", alpha=.2, color="red") +
  geom_rug(data=filter(dat_1_5, age=="Children" & sev.or.inte==1 & serotype=="DENV-1"), sides="t", alpha=.2, color="red") +
  geom_rug(data=filter(dat_1_5, age=="Adults" & sev.or.inte==0 & serotype=="Others"), sides="b", alpha=.2, color="blue") +
  geom_rug(data=filter(dat_1_5, age=="Adults" & sev.or.inte==1 & serotype=="Others"), sides="t", alpha=.2, color="blue") +
  scale_color_manual(values=c("red","blue"), labels=c("DENV-1","Others")) +
  scale_fill_manual(values=c("red","blue"), labels=c("DENV-1","Others")) +
  ylab("Odds ratio") +
  theme(legend.position="top", legend.title=element_blank(), axis.title.x=element_blank(), axis.text.x=element_text(angle=45)) +
  facet_grid_sc(cols=vars(biomarker), rows=vars(age), 
                scales=list(x=scales_x, y=scales_y1),
                labeller = as_labeller(c(age, biomarkers)))

# Plot for the last 5 biomarkers
p.6.10 <- ggplot(dat_6_10, aes(x=value)) +
  geom_ribbon(aes(ymin=lower, ymax=upper, fill=serotype), alpha=.15) +
  geom_vline(aes(xintercept=vline), color="black", linetype="dashed", alpha=.4) +
  geom_hline(yintercept=0, color="black", linetype="dashed", alpha=.4) +
  geom_line(aes(y=yhat, color=serotype), size=.5, alpha=.7) +
  geom_rug(data=filter(dat_6_10, age=="Children" & sev.or.inte==0 & serotype=="DENV-1"), sides="b", alpha=.2, color="red") +
  geom_rug(data=filter(dat_6_10, age=="Children" & sev.or.inte==1 & serotype=="DENV-1"), sides="t", alpha=.2, color="red") +
  geom_rug(data=filter(dat_6_10, age=="Adults" & sev.or.inte==0 & serotype=="Others"), sides="b", alpha=.2, color="blue") +
  geom_rug(data=filter(dat_6_10, age=="Adults" & sev.or.inte==1 & serotype=="Others"), sides="t", alpha=.2, color="blue") +
  scale_color_manual(values=c("red","blue"), labels=c("DENV-1","Others")) +
  scale_fill_manual(values=c("red","blue"), labels=c("DENV-1","Others")) +
  ylab("Odds ratio") +
  theme(legend.position="none", axis.title.x=element_blank(), axis.text.x=element_text(angle=45)) +
  facet_grid_sc(cols=vars(biomarker), rows=vars(age), 
                scales=list(x=scales_x, y=scales_y2),
                labeller = as_labeller(c(age, biomarkers)))

# Merge plots
p51 <- gridExtra::grid.arrange(p.1.5, p.6.10, nrow=2, heights=c(1,.86))
```

![](Elife-ERA-codes_files/figure-gfm/fig%20A5.1-1.png)<!-- -->

``` r
#ggsave(filename="A5_fig1.pdf", p3, dpi=300, width=8.5, height=6.2)
```

########################################################################################### 

#### Appendix 5—figure 2. Results from global model for severe/moderate dengue with the interaction with serotype.

``` r
# Set datadist for 'lrm' function [rms]
dd <- datadist(dat1); options(datadist="dd")

# For DENV-1
dd$limits["Adjust to","Serotype"] <- "DENV-1"

## For children
dd$limits["Adjust to","Age"] <- 10
dat_m1 <- get_pred2s(out="sev.or.inte", bio="VCAM", age=10, serotype="DENV-1", dat=dat1) %>% rename(value=VCAM) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='DENV-1')) %>% mutate(age="Children")
dat_m2 <- get_pred2s(out="sev.or.inte", bio="SDC", age=10, serotype="DENV-1", dat=dat1) %>% rename(value=SDC) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='DENV-1')) %>% mutate(age="Children")
dat_m3 <- get_pred2s(out="sev.or.inte", bio="Ang", age=10, serotype="DENV-1", dat=dat1) %>% rename(value=Ang) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='DENV-1')) %>% mutate(age="Children")
dat_m4 <- get_pred2s(out="sev.or.inte", bio="IL8", age=10, serotype="DENV-1", dat=dat1) %>% rename(value=IL8) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='DENV-1')) %>% mutate(age="Children")
dat_m5 <- get_pred2s(out="sev.or.inte", bio="IP10", age=10, serotype="DENV-1", dat=dat1) %>% rename(value=IP10) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='DENV-1')) %>% mutate(age="Children")
dat_m6 <- get_pred2s(out="sev.or.inte", bio="IL1RA", age=10, serotype="DENV-1", dat=dat1) %>% rename(value=IL1RA) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='DENV-1')) %>% mutate(age="Children")
dat_m7 <- get_pred2s(out="sev.or.inte", bio="CD163", age=10, serotype="DENV-1", dat=dat1) %>% rename(value=CD163) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='DENV-1')) %>% mutate(age="Children")
dat_m8 <- get_pred2s(out="sev.or.inte", bio="TREM", age=10, serotype="DENV-1", dat=dat1) %>% rename(value=TREM) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='DENV-1')) %>% mutate(age="Children")
dat_m9 <- get_pred2s(out="sev.or.inte", bio="Fer", age=10, serotype="DENV-1", dat=dat1) %>% rename(value=Fer) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='DENV-1')) %>% mutate(age="Children")
dat_m10 <- get_pred2s(out="sev.or.inte", bio="CRP", age=10, serotype="DENV-1", dat=dat1) %>% rename(value=CRP) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='DENV-1')) %>% mutate(age="Children")

## For adults
dd$limits["Adjust to","Age"] <- 25
dat_mg1 <- get_pred2s(out="sev.or.inte", bio="VCAM", age=25, serotype="DENV-1", dat=dat1) %>% rename(value=VCAM) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='DENV-1')) %>% mutate(age="Adults")
dat_mg2 <- get_pred2s(out="sev.or.inte", bio="SDC", age=25, serotype="DENV-1", dat=dat1) %>% rename(value=SDC) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='DENV-1')) %>% mutate(age="Adults")
dat_mg3 <- get_pred2s(out="sev.or.inte", bio="Ang", age=25, serotype="DENV-1", dat=dat1) %>% rename(value=Ang) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='DENV-1')) %>% mutate(age="Adults")
dat_mg4 <- get_pred2s(out="sev.or.inte", bio="IL8", age=25, serotype="DENV-1", dat=dat1) %>% rename(value=IL8) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='DENV-1')) %>% mutate(age="Adults")
dat_mg5 <- get_pred2s(out="sev.or.inte", bio="IP10", age=25, serotype="DENV-1", dat=dat1) %>% rename(value=IP10) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='DENV-1')) %>% mutate(age="Adults")
dat_mg6 <- get_pred2s(out="sev.or.inte", bio="IL1RA", age=25, serotype="DENV-1", dat=dat1) %>% rename(value=IL1RA) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='DENV-1')) %>% mutate(age="Adults")
dat_mg7 <- get_pred2s(out="sev.or.inte", bio="CD163", age=25, serotype="DENV-1", dat=dat1) %>% rename(value=CD163) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='DENV-1')) %>% mutate(age="Adults")
dat_mg8 <- get_pred2s(out="sev.or.inte", bio="TREM", age=25, serotype="DENV-1", dat=dat1) %>% rename(value=TREM) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='DENV-1')) %>% mutate(age="Adults")
dat_mg9 <- get_pred2s(out="sev.or.inte", bio="Fer", age=25, serotype="DENV-1", dat=dat1) %>% rename(value=Fer) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='DENV-1')) %>% mutate(age="Adults")
dat_mg10 <- get_pred2s(out="sev.or.inte", bio="CRP", age=25, serotype="DENV-1", dat=dat1) %>% rename(value=CRP) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='DENV-1')) %>% mutate(age="Adults")

dat_p1 <- rbind(dat_m1,dat_mg1, dat_m2,dat_mg2, dat_m3,dat_mg3, dat_m4,dat_mg4, dat_m5,dat_mg5, dat_m6,dat_mg6, dat_m7,dat_mg7, dat_m8,dat_mg8, dat_m9,dat_mg9, dat_m10,dat_mg10) %>%
  mutate(serotype = "DENV-1")

# For others
dd$limits["Adjust to","Serotype"] <- "Others"

## For children
dd$limits["Adjust to","Age"] <- 10
dat_m1 <- get_pred2s(out="sev.or.inte", bio="VCAM", age=10, serotype="Others", dat=dat1) %>% rename(value=VCAM) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='Others')) %>% mutate(age="Children")
dat_m2 <- get_pred2s(out="sev.or.inte", bio="SDC", age=10, serotype="Others", dat=dat1) %>% rename(value=SDC) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='Others')) %>% mutate(age="Children")
dat_m3 <- get_pred2s(out="sev.or.inte", bio="Ang", age=10, serotype="Others", dat=dat1) %>% rename(value=Ang) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='Others')) %>% mutate(age="Children")
dat_m4 <- get_pred2s(out="sev.or.inte", bio="IL8", age=10, serotype="Others", dat=dat1) %>% rename(value=IL8) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='Others')) %>% mutate(age="Children")
dat_m5 <- get_pred2s(out="sev.or.inte", bio="IP10", age=10, serotype="Others", dat=dat1) %>% rename(value=IP10) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='Others')) %>% mutate(age="Children")
dat_m6 <- get_pred2s(out="sev.or.inte", bio="IL1RA", age=10, serotype="Others", dat=dat1) %>% rename(value=IL1RA) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='Others')) %>% mutate(age="Children")
dat_m7 <- get_pred2s(out="sev.or.inte", bio="CD163", age=10, serotype="Others", dat=dat1) %>% rename(value=CD163) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='Others')) %>% mutate(age="Children")
dat_m8 <- get_pred2s(out="sev.or.inte", bio="TREM", age=10, serotype="Others", dat=dat1) %>% rename(value=TREM) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='Others')) %>% mutate(age="Children")
dat_m9 <- get_pred2s(out="sev.or.inte", bio="Fer", age=10, serotype="Others", dat=dat1) %>% rename(value=Fer) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='Others')) %>% mutate(age="Children")
dat_m10 <- get_pred2s(out="sev.or.inte", bio="CRP", age=10, serotype="Others", dat=dat1) %>% rename(value=CRP) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='Others')) %>% mutate(age="Children")

## For adults
dd$limits["Adjust to","Age"] <- 25
dat_mg1 <- get_pred2s(out="sev.or.inte", bio="VCAM", age=25, serotype="Others", dat=dat1) %>% rename(value=VCAM) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='Others')) %>% mutate(age="Adults")
dat_mg2 <- get_pred2s(out="sev.or.inte", bio="SDC", age=25, serotype="Others", dat=dat1) %>% rename(value=SDC) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='Others')) %>% mutate(age="Adults")
dat_mg3 <- get_pred2s(out="sev.or.inte", bio="Ang", age=25, serotype="Others", dat=dat1) %>% rename(value=Ang) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='Others')) %>% mutate(age="Adults")
dat_mg4 <- get_pred2s(out="sev.or.inte", bio="IL8", age=25, serotype="Others", dat=dat1) %>% rename(value=IL8) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='Others')) %>% mutate(age="Adults")
dat_mg5 <- get_pred2s(out="sev.or.inte", bio="IP10", age=25, serotype="Others", dat=dat1) %>% rename(value=IP10) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='Others')) %>% mutate(age="Adults")
dat_mg6 <- get_pred2s(out="sev.or.inte", bio="IL1RA", age=25, serotype="Others", dat=dat1) %>% rename(value=IL1RA) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='Others')) %>% mutate(age="Adults")
dat_mg7 <- get_pred2s(out="sev.or.inte", bio="CD163", age=25, serotype="Others", dat=dat1) %>% rename(value=CD163) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='Others')) %>% mutate(age="Adults")
dat_mg8 <- get_pred2s(out="sev.or.inte", bio="TREM", age=25, serotype="Others", dat=dat1) %>% rename(value=TREM) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='Others')) %>% mutate(age="Adults")
dat_mg9 <- get_pred2s(out="sev.or.inte", bio="Fer", age=25, serotype="Others", dat=dat1) %>% rename(value=Fer) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='Others')) %>% mutate(age="Adults")
dat_mg10 <- get_pred2s(out="sev.or.inte", bio="CRP", age=25, serotype="Others", dat=dat1) %>% rename(value=CRP) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$Serotype=='Others')) %>% mutate(age="Adults")

dat_p2 <- rbind(dat_m1,dat_mg1, dat_m2,dat_mg2, dat_m3,dat_mg3, dat_m4,dat_mg4, dat_m5,dat_mg5, dat_m6,dat_mg6, dat_m7,dat_mg7, dat_m8,dat_mg8, dat_m9,dat_mg9, dat_m10,dat_mg10) %>%
  mutate(serotype = "Others")

# Merge data for children and adults
vis1 <- rbind(dat_p1, dat_p2) %>% mutate(outcome = "outcome 1")

# Merge data for plots
tmp0 <- dat1 %>% filter(!is.na(Serotype)) %>% arrange(Serotype) %>% select(sev.or.inte)
tmp <- data.frame(sev.or.inte = rep(tmp0$sev.or.inte, 20))
  
tmp_p1 <- vis1 %>%
  arrange(biomarker, age) %>%
  mutate(value1 = 2^value,
         serotype = factor(serotype, levels=c("DENV-1", "Others")),
         age = factor(age, levels=c("Children", "Adults"))) %>%
  bind_cols(., tmp)

vline <- tmp_p1 %>%
  filter(yhat==0) %>%
  select(biomarker, value, value1) %>%
  rename(vline=value, vline1=value1)

dat_p <- left_join(tmp_p1, vline, by="biomarker")

# Alternative data to limit to 5th - 95th of the values
dat_alt <- dat_p %>%
  group_by(biomarker, model) %>%
  mutate(up = quantile(value, .95),
         lo = quantile(value, .05)) %>%
  ungroup() %>%
  mutate(is.outlier = value<lo | value>up | value<0,
         value = ifelse(is.outlier, NA, value),
         value1 = ifelse(is.outlier, NA, value1),
         yhat = ifelse(is.outlier, NA, yhat),
         lower = ifelse(is.outlier, NA, lower),
         upper = ifelse(is.outlier, NA, upper))

dat_1_5 <- dat_alt %>% 
  filter(biomarker %in% c("VCAM", "SDC", "Ang", "IL8", "IP10")) %>%
  mutate(biomarker = factor(biomarker, levels=c("VCAM", "SDC", "Ang", "IL8", "IP10")))

dat_6_10 <- dat_alt %>% 
  filter(biomarker %in% c("IL1RA", "CD163", "TREM", "Fer", "CRP")) %>%
  mutate(biomarker = factor(biomarker, levels=c("IL1RA", "CD163", "TREM", "Fer", "CRP")))

# Modify facets' scales
#devtools::install_github("zeehio/facetscales")
library(facetscales)

xVCAM <- c(1,4,15,60,250,1000,4000)
xSDC <- c(1400,2000,2800,4000,5600)
xAng <- c(50,100,250,500,1000,2000)
xIL8 <- c(5,7,10,14,20,28,40)
xIP10 <- c(25,100,400,1600,6400)
xIL1RA <- c(1000,2000,4000,8000,16000)
xCD163 <- c(75,150,300,600)
xTREM <- c(35,50,70,100,140,200)
xFer <- c(50,100,200,400,800)
xCRP <- c(2.5,5,10,20,40,80)

scales_x <- list(
  `VCAM` = scale_x_continuous(breaks = log2(xVCAM), labels = xVCAM),
  `SDC` = scale_x_continuous(breaks = log2(xSDC), labels = xSDC),
  `Ang` = scale_x_continuous(breaks = log2(xAng), labels = xAng), 
  `IL8` = scale_x_continuous(breaks = log2(xIL8), labels = xIL8), 
  `IP10` = scale_x_continuous(breaks = log2(xIP10), labels = xIP10), 
  `IL1RA` = scale_x_continuous(breaks = log2(xIL1RA), labels = xIL1RA), 
  `CD163` = scale_x_continuous(breaks = log2(xCD163), labels = xCD163), 
  `TREM` = scale_x_continuous(breaks = log2(xTREM), labels = xTREM), 
  `Fer` = scale_x_continuous(breaks = log2(xFer), labels = xFer), 
  `CRP` = scale_x_continuous(breaks = log2(xCRP), labels = xCRP)
)

ybreak1 <- c(.125, .25, .5, 1, 2, 4, 8)
ybreak2 <- c(.125, .25, .5, 1, 2, 4)

scales_y1 <- list(
  `Children` = scale_y_continuous(limits = c(log(.124), log(8.01)), breaks = log(ybreak1), labels = ybreak1),
  `Adults` = scale_y_continuous(limits = c(log(.124), log(8.01)), breaks = log(ybreak1), labels = ybreak1)
)

scales_y2 <- list(
  `Children` = scale_y_continuous(limits = c(log(.124), log(4.01)), breaks = log(ybreak2), labels = ybreak2),
  `Adults` = scale_y_continuous(limits = c(log(.124), log(4.01)), breaks = log(ybreak2), labels = ybreak2)
)

# Set facets' names
age <- c(`Children` = "Children",
         `Adults` = "Adults")

biomarkers <- c(`VCAM` = "VCAM-1 (ng/ml)",
                `SDC`= "SDC-1 (pg/ml)",
                `Ang` = "Ang-2 (pg/ml)",
                `IL8` = "IL-8 (pg/ml)",
                `IP10` = "IP-10 (pg/ml)",
                `IL1RA` = "IL-1RA (pg/ml)",
                `CD163` = "sCD163 (ng/ml)",
                `TREM` = "sTREM-1 (pg/ml)",
                `Fer` = "Ferritin (ng/ml)",
                `CRP` = "CRP (mg/l)")

# Plot for the first 5 biomarkers
p.1.5 <- ggplot(dat_1_5, aes(x=value)) +
  geom_ribbon(aes(ymin=lower, ymax=upper, fill=serotype), alpha=.15) +
  geom_vline(aes(xintercept=vline), color="black", linetype="dashed", alpha=.4) +
  geom_hline(yintercept=0, color="black", linetype="dashed", alpha=.4) +
  geom_line(aes(y=yhat, color=serotype), size=.5, alpha=.7) +
  geom_rug(data=filter(dat_1_5, age=="Children" & sev.or.inte==0 & serotype=="DENV-1"), sides="b", alpha=.2, color="red") +
  geom_rug(data=filter(dat_1_5, age=="Children" & sev.or.inte==1 & serotype=="DENV-1"), sides="t", alpha=.2, color="red") +
  geom_rug(data=filter(dat_1_5, age=="Adults" & sev.or.inte==0 & serotype=="Others"), sides="b", alpha=.2, color="blue") +
  geom_rug(data=filter(dat_1_5, age=="Adults" & sev.or.inte==1 & serotype=="Others"), sides="t", alpha=.2, color="blue") +
  scale_color_manual(values=c("red","blue"), labels=c("DENV-1","Others")) +
  scale_fill_manual(values=c("red","blue"), labels=c("DENV-1","Others")) +
  ylab("Odds ratio") +
  theme(legend.position="top", legend.title=element_blank(), axis.title.x=element_blank(), axis.text.x=element_text(angle=45)) +
  facet_grid_sc(cols=vars(biomarker), rows=vars(age), 
                scales=list(x=scales_x, y=scales_y1),
                labeller = as_labeller(c(age, biomarkers)))

# Plot for the last 5 biomarkers
p.6.10 <- ggplot(dat_6_10, aes(x=value)) +
  geom_ribbon(aes(ymin=lower, ymax=upper, fill=serotype), alpha=.15) +
  geom_vline(aes(xintercept=vline), color="black", linetype="dashed", alpha=.4) +
  geom_hline(yintercept=0, color="black", linetype="dashed", alpha=.4) +
  geom_line(aes(y=yhat, color=serotype), size=.5, alpha=.7) +
  geom_rug(data=filter(dat_6_10, age=="Children" & sev.or.inte==0 & serotype=="DENV-1"), sides="b", alpha=.2, color="red") +
  geom_rug(data=filter(dat_6_10, age=="Children" & sev.or.inte==1 & serotype=="DENV-1"), sides="t", alpha=.2, color="red") +
  geom_rug(data=filter(dat_6_10, age=="Adults" & sev.or.inte==0 & serotype=="Others"), sides="b", alpha=.2, color="blue") +
  geom_rug(data=filter(dat_6_10, age=="Adults" & sev.or.inte==1 & serotype=="Others"), sides="t", alpha=.2, color="blue") +
  scale_color_manual(values=c("red","blue"), labels=c("DENV-1","Others")) +
  scale_fill_manual(values=c("red","blue"), labels=c("DENV-1","Others")) +
  ylab("Odds ratio") +
  theme(legend.position="none", axis.title.x=element_blank(), axis.text.x=element_text(angle=45)) +
  facet_grid_sc(cols=vars(biomarker), rows=vars(age), 
                scales=list(x=scales_x, y=scales_y2),
                labeller = as_labeller(c(age, biomarkers)))

# Merge plots
p52 <- gridExtra::grid.arrange(p.1.5, p.6.10, nrow=2, heights=c(1,.86))
```

![](Elife-ERA-codes_files/figure-gfm/fig%20A5.2-1.png)<!-- -->

``` r
#ggsave(filename="A5_fig2.pdf", p52, dpi=300, width=8.5, height=6.2)
```

########################################################################################### 

#### Appendix 5—table 1. Results from single models for severe/moderate dengue with the interaction with serotype.

``` r
# Set datadist for 'lrm' function [rms]
dd <- datadist(dat1); options(datadist="dd")

# Get ORs and 95% CIs
tmp1 <- data.frame(
  sort1 = rep(c(1:10), 2),
  sort2 = c(rep(2,10), rep(3,10)),
  bio = rep(c("VCAM", "SDC", "Ang", "IL8", "IP10", "IL1RA", "CD163", "TREM", "Fer", "CRP"), 2),
  ref1 = c(ref0-1, ref0),
  ref2 = c(ref0, ref0+1)
) %>%
  arrange(sort1, sort2) %>%
  group_by(sort1, sort2) %>%
  do(cbind(.,
           # DENV-1 - children
           or1c = get_est1s(out="sev.or.inte", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="OR", age=10, serotype=1, dat=dat1),
           lo1c = get_est1s(out="sev.or.inte", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="loCI", age=10, serotype=1, dat=dat1),
           up1c = get_est1s(out="sev.or.inte", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="upCI", age=10, serotype=1, dat=dat1),
           # Others - children
           or1a = get_est1s(out="sev.or.inte", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="OR", age=10, serotype=2, dat=dat1),
           lo1a = get_est1s(out="sev.or.inte", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="loCI", age=10, serotype=2, dat=dat1),
           up1a = get_est1s(out="sev.or.inte", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="upCI", age=10, serotype=2, dat=dat1),
           # DENV-1 - adults
           or2c = get_est1s(out="sev.or.inte", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="OR", age=25, serotype=1, dat=dat1),
           lo2c = get_est1s(out="sev.or.inte", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="loCI", age=25, serotype=1, dat=dat1),
           up2c = get_est1s(out="sev.or.inte", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="upCI", age=25, serotype=1, dat=dat1),
           # Others - adults
           or2a = get_est1s(out="sev.or.inte", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="OR", age=25, serotype=2, dat=dat1),
           lo2a = get_est1s(out="sev.or.inte", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="loCI", age=25, serotype=2, dat=dat1),
           up2a = get_est1s(out="sev.or.inte", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="upCI", age=25, serotype=2, dat=dat1))) %>%
  ungroup()
for (i in 6:17) {tmp1[[i]] <- sprintf("%.2f", round(tmp1[[i]],2))}

# Get p-values
tmp2 <- data.frame(
  sort1 = c(1:10),
  sort2 = rep(1,10),
  bio = c("VCAM", "SDC", "Ang", "IL8", "IP10", "IL1RA", "CD163", "TREM", "Fer", "CRP")
) %>%
  group_by(sort1) %>%
  do(cbind(.,
           p.s1 = get_est1s(out="sev.or.inte", bio=.$bio, est="p", dat=dat1),
           p.s1_int = get_est1s(out="sev.or.inte", bio=.$bio, est="p int", dat=dat1),
           p.s1_int_age = get_est1s(out="sev.or.inte", bio=.$bio, est="p int age", dat=dat1),
           p.s1_int_serotype = get_est1s(out="sev.or.inte", bio=.$bio, est="p int serotype", dat=dat1))) %>%
  ungroup()
for (i in 4:7) {tmp2[[i]] <- ifelse(tmp2[[i]]<0.001, "<0.001", sprintf("%.3f", round(tmp2[[i]],3)))}  

# Combine results
res1 <- bind_rows(tmp1, tmp2) %>%
  arrange(sort1, sort2) %>%
  mutate(bio = ifelse(is.na(ref1), as.character(bio),
                      ifelse(bio!="Vir", paste(" -", round(2^ref2,0), "vs", round(2^ref1,0), sep=" "), paste(" -", round(10^ref2,0), "vs", round(10^ref1,0), sep=" "))),
         or.sc1 = ifelse(is.na(lo1c), NA, paste(or1c, " (", lo1c, "-", up1c, ")", sep="")),
         or.sa1 = ifelse(is.na(lo1a), NA, paste(or1a, " (", lo1a, "-", up1a, ")", sep="")),
         or.sc2 = ifelse(is.na(lo2c), NA, paste(or2c, " (", lo2c, "-", up2c, ")", sep="")),
         or.sa2 = ifelse(is.na(lo2a), NA, paste(or2a, " (", lo2a, "-", up2a, ")", sep=""))) %>%
  select(bio, or.sc1, or.sa1, or.sc2, or.sa2, p.s1, p.s1_int, p.s1_int_age, p.s1_int_serotype)

names(res1) <- c("", "OR (children - DENV-1)", "OR (children - others)", "OR (adults - DENV-1)", "OR (adults - others)",
                 "P overall", "P overall interaction", "P interaction with age", "P interaction with serotype")

# Report results
knitr::kable(res1)
```

|                  | OR (children - DENV-1) | OR (children - others) | OR (adults - DENV-1) | OR (adults - others) | P overall | P overall interaction | P interaction with age | P interaction with serotype |
|:-----------------|:-----------------------|:-----------------------|:---------------------|:---------------------|:----------|:----------------------|:-----------------------|:----------------------------|
| VCAM             |                        |                        |                      |                      | 0.008     | 0.822                 | 0.729                  | 0.565                       |
| \- 1636 vs 818   | 1.22 (1.02-1.46)       | 1.14 (0.95-1.37)       | 1.42 (1.13-1.77)     | 1.32 (1.10-1.58)     |           |                       |                        |                             |
| \- 3272 vs 1636  | 1.28 (0.99-1.66)       | 1.18 (0.91-1.54)       | 1.57 (1.15-2.14)     | 1.45 (1.13-1.86)     |           |                       |                        |                             |
| SDC              |                        |                        |                      |                      | &lt;0.001 | 0.087                 | 0.326                  | 0.352                       |
| \- 2519 vs 1260  | 2.35 (0.89-6.25)       | 2.56 (0.97-6.74)       | 2.85 (0.80-10.12)    | 3.10 (1.11-8.64)     |           |                       |                        |                             |
| \- 5039 vs 2519  | 2.41 (1.47-3.96)       | 1.22 (0.69-2.13)       | 6.39 (2.81-14.53)    | 3.23 (1.73-6.01)     |           |                       |                        |                             |
| Ang              |                        |                        |                      |                      | &lt;0.001 | 0.935                 | 0.923                  | 0.702                       |
| \- 1204 vs 602   | 1.67 (1.35-2.05)       | 1.51 (1.21-1.88)       | 1.64 (1.28-2.10)     | 1.48 (1.20-1.83)     |           |                       |                        |                             |
| \- 2409 vs 1204  | 2.46 (1.60-3.79)       | 1.97 (1.25-3.11)       | 2.40 (1.44-4.02)     | 1.92 (1.28-2.89)     |           |                       |                        |                             |
| IL8              |                        |                        |                      |                      | &lt;0.001 | &lt;0.001             | &lt;0.001              | 0.104                       |
| \- 14 vs 7       | 1.52 (1.00-2.31)       | 1.12 (0.76-1.65)       | 2.83 (1.57-5.12)     | 2.09 (1.34-3.26)     |           |                       |                        |                             |
| \- 28 vs 14      | 1.22 (0.90-1.65)       | 0.84 (0.60-1.17)       | 3.04 (1.88-4.94)     | 2.10 (1.45-3.04)     |           |                       |                        |                             |
| IP10             |                        |                        |                      |                      | &lt;0.001 | 0.975                 | 0.950                  | 0.681                       |
| \- 3093 vs 1546  | 1.54 (1.28-1.84)       | 1.39 (1.13-1.70)       | 1.64 (1.28-2.10)     | 1.48 (1.20-1.83)     |           |                       |                        |                             |
| \- 6186 vs 3093  | 1.84 (1.39-2.43)       | 1.60 (1.17-2.20)       | 2.00 (1.37-2.92)     | 1.75 (1.27-2.40)     |           |                       |                        |                             |
| IL1RA            |                        |                        |                      |                      | &lt;0.001 | 0.577                 | 0.280                  | 0.805                       |
| \- 6434 vs 3217  | 1.68 (1.29-2.18)       | 1.51 (1.20-1.90)       | 1.96 (1.34-2.86)     | 1.76 (1.29-2.38)     |           |                       |                        |                             |
| \- 12868 vs 6434 | 1.87 (1.43-2.44)       | 1.77 (1.28-2.46)       | 1.73 (1.23-2.44)     | 1.64 (1.18-2.29)     |           |                       |                        |                             |
| CD163            |                        |                        |                      |                      | 0.002     | 0.983                 | 0.831                  | 0.932                       |
| \- 295 vs 147    | 1.50 (1.03-2.20)       | 1.60 (1.03-2.48)       | 1.59 (1.02-2.47)     | 1.69 (1.18-2.42)     |           |                       |                        |                             |
| \- 589 vs 295    | 1.35 (0.94-1.94)       | 1.48 (0.99-2.21)       | 1.50 (0.92-2.43)     | 1.64 (1.04-2.59)     |           |                       |                        |                             |
| TREM             |                        |                        |                      |                      | 0.146     | 0.979                 | 0.998                  | 0.597                       |
| \- 85 vs 42      | 1.55 (0.94-2.57)       | 2.03 (1.22-3.39)       | 1.68 (0.92-3.07)     | 2.20 (1.21-4.02)     |           |                       |                        |                             |
| \- 169 vs 85     | 1.08 (0.80-1.45)       | 1.16 (0.87-1.54)       | 1.08 (0.68-1.71)     | 1.16 (0.83-1.60)     |           |                       |                        |                             |
| Fer              |                        |                        |                      |                      | 0.112     | 0.139                 | 0.177                  | 0.711                       |
| \- 243 vs 122    | 1.18 (0.97-1.44)       | 1.12 (0.91-1.38)       | 1.14 (0.88-1.49)     | 1.08 (0.88-1.32)     |           |                       |                        |                             |
| \- 487 vs 243    | 1.28 (0.96-1.71)       | 1.32 (0.91-1.92)       | 0.87 (0.55-1.38)     | 0.90 (0.64-1.27)     |           |                       |                        |                             |
| CRP              |                        |                        |                      |                      | &lt;0.001 | 0.080                 | 0.029                  | 0.755                       |
| \- 28 vs 14      | 1.23 (1.06-1.43)       | 1.36 (1.12-1.66)       | 1.21 (0.94-1.57)     | 1.34 (1.07-1.69)     |           |                       |                        |                             |
| \- 56 vs 28      | 1.06 (0.85-1.33)       | 1.24 (0.97-1.58)       | 1.25 (0.93-1.68)     | 1.46 (1.13-1.87)     |           |                       |                        |                             |

########################################################################################### 

#### Appendix 5—table 2. Results from global model for severe/moderate dengue with the interaction with serotype.

``` r
# Set datadist for 'lrm' function [rms]
dd <- datadist(dat1); options(datadist="dd")

# Get ORs and 95% CIs
tmp1 <- data.frame(
  sort1 = rep(c(1:10), 2),
  sort2 = c(rep(2,10), rep(3,10)),
  bio = rep(c("VCAM", "SDC", "Ang", "IL8", "IP10", "IL1RA", "CD163", "TREM", "Fer", "CRP"), 2),
  ref1 = c(ref0-1, ref0),
  ref2 = c(ref0, ref0+1)
) %>%
  arrange(sort1, sort2) %>%
  group_by(sort1, sort2) %>%
  do(cbind(.,
           # DENV-1 - children
           or1c = get_est2s(out="sev.or.inte", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="OR", age=10, serotype=1, dat=dat1),
           lo1c = get_est2s(out="sev.or.inte", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="loCI", age=10, serotype=1, dat=dat1),
           up1c = get_est2s(out="sev.or.inte", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="upCI", age=10, serotype=1, dat=dat1),
           # Others - children
           or1a = get_est2s(out="sev.or.inte", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="OR", age=10, serotype=2, dat=dat1),
           lo1a = get_est2s(out="sev.or.inte", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="loCI", age=10, serotype=2, dat=dat1),
           up1a = get_est2s(out="sev.or.inte", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="upCI", age=10, serotype=2, dat=dat1),
           # DENV-1 - adults
           or2c = get_est2s(out="sev.or.inte", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="OR", age=25, serotype=1, dat=dat1),
           lo2c = get_est2s(out="sev.or.inte", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="loCI", age=25, serotype=1, dat=dat1),
           up2c = get_est2s(out="sev.or.inte", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="upCI", age=25, serotype=1, dat=dat1),
           # Others - adults
           or2a = get_est2s(out="sev.or.inte", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="OR", age=25, serotype=2, dat=dat1),
           lo2a = get_est2s(out="sev.or.inte", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="loCI", age=25, serotype=2, dat=dat1),
           up2a = get_est2s(out="sev.or.inte", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="upCI", age=25, serotype=2, dat=dat1))) %>%
  ungroup()
for (i in 6:17) {tmp1[[i]] <- sprintf("%.2f", round(tmp1[[i]],2))}

# Get p-values
tmp2 <- data.frame(
  sort1 = c(1:10),
  sort2 = rep(1,10),
  bio = c("VCAM", "SDC", "Ang", "IL8", "IP10", "IL1RA", "CD163", "TREM", "Fer", "CRP")
) %>%
  group_by(sort1) %>%
  do(cbind(.,
           p.s2 = get_est2s(out="sev.or.inte", bio=.$bio, est="p", dat=dat1),
           p.s2_int = get_est2s(out="sev.or.inte", bio=.$bio, est="p int", dat=dat1),
           p.s2_int_age = get_est2s(out="sev.or.inte", bio=.$bio, est="p int age", dat=dat1),
           p.s2_int_serotype = get_est2s(out="sev.or.inte", bio=.$bio, est="p int serotype", dat=dat1))) %>%
  ungroup()
for (i in 4:7) {tmp2[[i]] <- ifelse(tmp2[[i]]<0.001, "<0.001", sprintf("%.3f", round(tmp2[[i]],3)))}  

# Combine results
res2 <- bind_rows(tmp1, tmp2) %>%
  arrange(sort1, sort2) %>%
  mutate(bio = ifelse(is.na(ref1), as.character(bio),
                      ifelse(bio!="Vir", paste(" -", round(2^ref2,0), "vs", round(2^ref1,0), sep=" "), paste(" -", round(10^ref2,0), "vs", round(10^ref1,0), sep=" "))),
         or.sc1 = ifelse(is.na(lo1c), NA, paste(or1c, " (", lo1c, "-", up1c, ")", sep="")),
         or.sa1 = ifelse(is.na(lo1a), NA, paste(or1a, " (", lo1a, "-", up1a, ")", sep="")),
         or.sc2 = ifelse(is.na(lo2c), NA, paste(or2c, " (", lo2c, "-", up2c, ")", sep="")),
         or.sa2 = ifelse(is.na(lo2a), NA, paste(or2a, " (", lo2a, "-", up2a, ")", sep=""))) %>%
  select(bio, or.sc1, or.sa1, or.sc2, or.sa2, p.s2, p.s2_int, p.s2_int_age, p.s2_int_serotype)

names(res2) <- c("", "OR (children - DENV-1)", "OR (children - others)", "OR (adults - DENV-1)", "OR (adults - others)",
                 "P overall", "P overall interaction", "P interaction with age", "P interaction with serotype")

# Report results
knitr::kable(res2)
```

|                  | OR (children - DENV-1) | OR (children - others) | OR (adults - DENV-1) | OR (adults - others) | P overall | P overall interaction | P interaction with age | P interaction with serotype |
|:-----------------|:-----------------------|:-----------------------|:---------------------|:---------------------|:----------|:----------------------|:-----------------------|:----------------------------|
| VCAM             |                        |                        |                      |                      | 0.449     | 0.258                 | 0.248                  | 0.327                       |
| \- 1636 vs 818   | 0.84 (0.62-1.13)       | 0.88 (0.66-1.17)       | 1.18 (0.81-1.72)     | 1.24 (0.93-1.65)     |           |                       |                        |                             |
| \- 3272 vs 1636  | 0.75 (0.50-1.12)       | 0.85 (0.58-1.25)       | 1.14 (0.70-1.85)     | 1.29 (0.88-1.90)     |           |                       |                        |                             |
| SDC              |                        |                        |                      |                      | 0.027     | 0.788                 | 0.821                  | 0.316                       |
| \- 2519 vs 1260  | 3.21 (0.79-12.94)      | 1.38 (0.38-4.94)       | 5.98 (1.00-35.72)    | 2.57 (0.66-10.01)    |           |                       |                        |                             |
| \- 5039 vs 2519  | 2.82 (1.21-6.57)       | 1.45 (0.62-3.40)       | 3.70 (1.08-12.72)    | 1.90 (0.80-4.52)     |           |                       |                        |                             |
| Ang              |                        |                        |                      |                      | 0.067     | 0.102                 | 0.043                  | 0.472                       |
| \- 1204 vs 602   | 1.65 (1.10-2.47)       | 1.79 (1.18-2.72)       | 0.89 (0.55-1.44)     | 0.97 (0.67-1.40)     |           |                       |                        |                             |
| \- 2409 vs 1204  | 2.22 (1.21-4.05)       | 1.73 (0.94-3.17)       | 1.25 (0.63-2.48)     | 0.97 (0.58-1.62)     |           |                       |                        |                             |
| IL8              |                        |                        |                      |                      | &lt;0.001 | &lt;0.001             | &lt;0.001              | 0.591                       |
| \- 14 vs 7       | 0.86 (0.48-1.51)       | 0.73 (0.44-1.24)       | 2.14 (1.01-4.53)     | 1.84 (1.05-3.20)     |           |                       |                        |                             |
| \- 28 vs 14      | 0.68 (0.41-1.13)       | 0.49 (0.29-0.82)       | 2.60 (1.27-5.32)     | 1.88 (1.20-2.96)     |           |                       |                        |                             |
| IP10             |                        |                        |                      |                      | 0.068     | 0.875                 | 0.715                  | 0.888                       |
| \- 3093 vs 1546  | 0.98 (0.69-1.40)       | 0.91 (0.65-1.26)       | 0.86 (0.53-1.40)     | 0.80 (0.53-1.19)     |           |                       |                        |                             |
| \- 6186 vs 3093  | 1.27 (0.77-2.10)       | 1.10 (0.70-1.74)       | 1.04 (0.54-2.00)     | 0.90 (0.53-1.53)     |           |                       |                        |                             |
| IL1RA            |                        |                        |                      |                      | &lt;0.001 | 0.333                 | 0.230                  | 0.711                       |
| \- 6434 vs 3217  | 2.63 (1.65-4.20)       | 2.11 (1.32-3.39)       | 2.50 (1.28-4.88)     | 2.01 (1.15-3.51)     |           |                       |                        |                             |
| \- 12868 vs 6434 | 2.38 (1.45-3.91)       | 1.90 (1.19-3.04)       | 1.65 (0.88-3.09)     | 1.32 (0.78-2.23)     |           |                       |                        |                             |
| CD163            |                        |                        |                      |                      | 0.340     | 0.661                 | 0.455                  | 0.769                       |
| \- 295 vs 147    | 1.20 (0.65-2.20)       | 1.54 (0.83-2.87)       | 1.24 (0.68-2.27)     | 1.61 (0.92-2.80)     |           |                       |                        |                             |
| \- 589 vs 295    | 1.20 (0.77-1.87)       | 1.20 (0.75-1.94)       | 1.49 (0.84-2.67)     | 1.50 (0.85-2.61)     |           |                       |                        |                             |
| TREM             |                        |                        |                      |                      | 0.441     | 0.289                 | 0.306                  | 0.071                       |
| \- 85 vs 42      | 0.90 (0.48-1.67)       | 1.23 (0.66-2.28)       | 0.84 (0.36-1.94)     | 1.15 (0.54-2.42)     |           |                       |                        |                             |
| \- 169 vs 85     | 0.56 (0.33-0.96)       | 1.18 (0.79-1.77)       | 0.34 (0.15-0.76)     | 0.71 (0.43-1.19)     |           |                       |                        |                             |
| Fer              |                        |                        |                      |                      | 0.067     | 0.033                 | 0.013                  | 0.331                       |
| \- 243 vs 122    | 1.36 (1.00-1.85)       | 1.25 (0.92-1.69)       | 0.92 (0.62-1.36)     | 0.84 (0.63-1.13)     |           |                       |                        |                             |
| \- 487 vs 243    | 1.08 (0.72-1.63)       | 1.59 (0.97-2.61)       | 0.55 (0.29-1.05)     | 0.81 (0.51-1.29)     |           |                       |                        |                             |
| CRP              |                        |                        |                      |                      | 0.156     | 0.103                 | 0.241                  | 0.136                       |
| \- 28 vs 14      | 0.95 (0.78-1.16)       | 1.25 (0.98-1.60)       | 0.94 (0.67-1.32)     | 1.24 (0.90-1.70)     |           |                       |                        |                             |
| \- 56 vs 28      | 0.80 (0.60-1.07)       | 1.07 (0.80-1.44)       | 1.13 (0.75-1.70)     | 1.51 (1.09-2.09)     |           |                       |                        |                             |

########################################################################################### 

#### Appendix 6—figure 1. Results from models for severe dengue endpoint.

``` r
# Set datadist for 'lrm' function [rms]
dd <- datadist(dat1); options(datadist="dd")

# For children
dd$limits["Adjust to","Age"] <- 10

dat_m1 <- get_pred1(out="sev.only", bio="VCAM", age=10, dat=dat1) %>% rename(value=VCAM) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_m2 <- get_pred1(out="sev.only", bio="SDC", age=10, dat=dat1) %>% rename(value=SDC) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_m3 <- get_pred1(out="sev.only", bio="Ang", age=10, dat=dat1) %>% rename(value=Ang) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_m4 <- get_pred1(out="sev.only", bio="IL8", age=10, dat=dat1) %>% rename(value=IL8) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_m5 <- get_pred1(out="sev.only", bio="IP10", age=10, dat=dat1) %>% rename(value=IP10) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_m6 <- get_pred1(out="sev.only", bio="IL1RA", age=10, dat=dat1) %>% rename(value=IL1RA) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_m7 <- get_pred1(out="sev.only", bio="CD163", age=10, dat=dat1) %>% rename(value=CD163) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_m8 <- get_pred1(out="sev.only", bio="TREM", age=10, dat=dat1) %>% rename(value=TREM) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_m9 <- get_pred1(out="sev.only", bio="Fer", age=10, dat=dat1) %>% rename(value=Fer) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_m10 <- get_pred1(out="sev.only", bio="CRP", age=10, dat=dat1) %>% rename(value=CRP) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))

dat_mg1 <- get_pred2(out="sev.only", bio="VCAM", age=10, dat=dat1) %>% rename(value=VCAM) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_mg2 <- get_pred2(out="sev.only", bio="SDC", age=10, dat=dat1) %>% rename(value=SDC) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_mg3 <- get_pred2(out="sev.only", bio="Ang", age=10, dat=dat1) %>% rename(value=Ang) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_mg4 <- get_pred2(out="sev.only", bio="IL8", age=10, dat=dat1) %>% rename(value=IL8) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_mg5 <- get_pred2(out="sev.only", bio="IP10", age=10, dat=dat1) %>% rename(value=IP10) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_mg6 <- get_pred2(out="sev.only", bio="IL1RA", age=10, dat=dat1) %>% rename(value=IL1RA) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_mg7 <- get_pred2(out="sev.only", bio="CD163", age=10, dat=dat1) %>% rename(value=CD163) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_mg8 <- get_pred2(out="sev.only", bio="TREM", age=10, dat=dat1) %>% rename(value=TREM) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_mg9 <- get_pred2(out="sev.only", bio="Fer", age=10, dat=dat1) %>% rename(value=Fer) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_mg10 <- get_pred2(out="sev.only", bio="CRP", age=10, dat=dat1) %>% rename(value=CRP) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))

dat_p1 <- rbind(dat_m1,dat_mg1, dat_m2,dat_mg2, dat_m3,dat_mg3, dat_m4,dat_mg4, dat_m5,dat_mg5, dat_m6,dat_mg6, dat_m7,dat_mg7, dat_m8,dat_mg8, dat_m9,dat_mg9, dat_m10,dat_mg10) %>%
  mutate(age = "10 years")

# For adults
dd$limits["Adjust to","Age"] <- 25

dat_m1 <- get_pred1(out="sev.only", bio="VCAM", age=25, dat=dat1) %>% rename(value=VCAM) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_m2 <- get_pred1(out="sev.only", bio="SDC", age=25, dat=dat1) %>% rename(value=SDC) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_m3 <- get_pred1(out="sev.only", bio="Ang", age=25, dat=dat1) %>% rename(value=Ang) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_m4 <- get_pred1(out="sev.only", bio="IL8", age=25, dat=dat1) %>% rename(value=IL8) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_m5 <- get_pred1(out="sev.only", bio="IP10", age=25, dat=dat1) %>% rename(value=IP10) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_m6 <- get_pred1(out="sev.only", bio="IL1RA", age=25, dat=dat1) %>% rename(value=IL1RA) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_m7 <- get_pred1(out="sev.only", bio="CD163", age=25, dat=dat1) %>% rename(value=CD163) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_m8 <- get_pred1(out="sev.only", bio="TREM", age=25, dat=dat1) %>% rename(value=TREM) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_m9 <- get_pred1(out="sev.only", bio="Fer", age=25, dat=dat1) %>% rename(value=Fer) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_m10 <- get_pred1(out="sev.only", bio="CRP", age=25, dat=dat1) %>% rename(value=CRP) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))

dat_mg1 <- get_pred2(out="sev.only", bio="VCAM", age=25, dat=dat1) %>% rename(value=VCAM) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_mg2 <- get_pred2(out="sev.only", bio="SDC", age=25, dat=dat1) %>% rename(value=SDC) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_mg3 <- get_pred2(out="sev.only", bio="Ang", age=25, dat=dat1) %>% rename(value=Ang) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_mg4 <- get_pred2(out="sev.only", bio="IL8", age=25, dat=dat1) %>% rename(value=IL8) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_mg5 <- get_pred2(out="sev.only", bio="IP10", age=25, dat=dat1) %>% rename(value=IP10) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_mg6 <- get_pred2(out="sev.only", bio="IL1RA", age=25, dat=dat1) %>% rename(value=IL1RA) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_mg7 <- get_pred2(out="sev.only", bio="CD163", age=25, dat=dat1) %>% rename(value=CD163) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_mg8 <- get_pred2(out="sev.only", bio="TREM", age=25, dat=dat1) %>% rename(value=TREM) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_mg9 <- get_pred2(out="sev.only", bio="Fer", age=25, dat=dat1) %>% rename(value=Fer) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_mg10 <- get_pred2(out="sev.only", bio="CRP", age=25, dat=dat1) %>% rename(value=CRP) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))

dat_p2 <- rbind(dat_m1,dat_mg1, dat_m2,dat_mg2, dat_m3,dat_mg3, dat_m4,dat_mg4, dat_m5,dat_mg5, dat_m6,dat_mg6, dat_m7,dat_mg7, dat_m8,dat_mg8, dat_m9,dat_mg9, dat_m10,dat_mg10) %>%
  mutate(age = "25 years")

# Merge results for children and adults
vis2 <- rbind(dat_p1, dat_p2) %>%
  mutate(yhat = ifelse(yhat>log(32), NA, ifelse(yhat<log(1/32), NA, yhat)),
         lower = ifelse(lower>log(32), log(32), ifelse(lower<log(1/32), log(1/32), lower)),
         upper = ifelse(upper>log(32), log(32), ifelse(upper<log(1/32), log(1/32), upper)),
         outcome = "outcome 2")

# Merge data for plots
tmp0 <- dat %>% arrange(age15) %>% select(sev.only)
tmp <- data.frame(sev.only = rep(tmp0$sev.only, 20))

tmp_p1 <- vis2 %>%
  mutate(value1 = 2^value,
         age = factor(age, levels=c("10 years", "25 years")),
         gr = ifelse(model=="Single model" & age=="10 years", 1,
                     ifelse(model=="Single model" & age=="25 years", 2,
                            ifelse(model=="Global model" & age=="10 years", 3, 4))),
         gr = factor(gr, levels=c(1:4), labels=c("Single children", "Single adults", "Global children", "Global adults")),
         model = factor(model, levels=c("Single model", "Global model"))) %>%
  bind_cols(., tmp)

vline <- tmp_p1 %>%
  filter(model=="Single model") %>%
  filter(yhat==0) %>%
  select(biomarker, value, value1) %>%
  rename(vline=value, vline1=value1)

dat_p <- left_join(tmp_p1, vline, by="biomarker")

# Alternative data to limit to 5th - 95th of the values
dat_alt <- dat_p %>%
  group_by(biomarker, model) %>%
  mutate(up = quantile(value, .95),
         lo = quantile(value, .05)) %>%
  ungroup() %>%
  mutate(is.outlier = value<lo | value>up | value<0,
         value = ifelse(is.outlier, NA, value),
         value1 = ifelse(is.outlier, NA, value1),
         yhat = ifelse(is.outlier, NA, yhat),
         lower = ifelse(is.outlier, NA, lower),
         upper = ifelse(is.outlier, NA, upper))

dat_1_5 <- dat_alt %>% 
  filter(biomarker %in% c("VCAM", "SDC", "Ang", "IL8", "IP10")) %>%
  mutate(biomarker = factor(biomarker, levels=c("VCAM", "SDC", "Ang", "IL8", "IP10")))

dat_6_10 <- dat_alt %>% 
  filter(biomarker %in% c("IL1RA", "CD163", "TREM", "Fer", "CRP")) %>%
  mutate(biomarker = factor(biomarker, levels=c("IL1RA", "CD163", "TREM", "Fer", "CRP")))

# Modify facets' scales
#library(facetscales)

xVCAM <- c(1,4,15,60,250,1000,4000)
xSDC <- c(1400,2000,2800,4000,5600)
xAng <- c(50,100,250,500,1000,2000)
xIL8 <- c(5,7,10,14,20,28,40)
xIP10 <- c(25,100,400,1600,6400)
xIL1RA <- c(1000,2000,4000,8000,16000)
xCD163 <- c(75,150,300,600)
xTREM <- c(35,50,70,100,140,200)
xFer <- c(50,100,200,400,800)
xCRP <- c(2.5,5,10,20,40,80)

scales_x <- list(
  `VCAM` = scale_x_continuous(breaks = log2(xVCAM), labels = xVCAM),
  `SDC` = scale_x_continuous(breaks = log2(xSDC), labels = xSDC),
  `Ang` = scale_x_continuous(breaks = log2(xAng), labels = xAng), 
  `IL8` = scale_x_continuous(breaks = log2(xIL8), labels = xIL8), 
  `IP10` = scale_x_continuous(breaks = log2(xIP10), labels = xIP10), 
  `IL1RA` = scale_x_continuous(breaks = log2(xIL1RA), labels = xIL1RA), 
  `CD163` = scale_x_continuous(breaks = log2(xCD163), labels = xCD163), 
  `TREM` = scale_x_continuous(breaks = log2(xTREM), labels = xTREM), 
  `Fer` = scale_x_continuous(breaks = log2(xFer), labels = xFer), 
  `CRP` = scale_x_continuous(breaks = log2(xCRP), labels = xCRP)
)

ybreak1 <- c(.0625,.125, .25, .5, 1, 2, 4, 8, 16)
ybreak2 <- c(.0625,.125, .25, .5, 1, 2, 4, 8, 16)

scales_y1 <- list(
  `Single model` = scale_y_continuous(limits = c(log(.0624), log(16.01)), breaks = log(ybreak1), labels = ybreak1),
  `Global model` = scale_y_continuous(limits = c(log(.0624), log(16.01)), breaks = log(ybreak1), labels = ybreak1)
)

scales_y2 <- list(
  `Single model` = scale_y_continuous(limits = c(log(.0624), log(16.01)), breaks = log(ybreak2), labels = ybreak2),
  `Global model` = scale_y_continuous(limits = c(log(.0624), log(16.01)), breaks = log(ybreak2), labels = ybreak2)
)

# Set facets' names
models <- c(`Single model` = "Single model",
            `Global model` = "Global model")

biomarkers <- c(`VCAM` = "VCAM-1 (ng/ml)",
                `SDC`= "SDC-1 (pg/ml)",
                `Ang` = "Ang-2 (pg/ml)",
                `IL8` = "IL-8 (pg/ml)",
                `IP10` = "IP-10 (pg/ml)",
                `IL1RA` = "IL-1RA (pg/ml)",
                `CD163` = "sCD163 (ng/ml)",
                `TREM` = "sTREM-1 (pg/ml)",
                `Fer` = "Ferritin (ng/ml)",
                `CRP` = "CRP (mg/l)")

# Plot for the first 5 biomarkers
p.1.5 <- ggplot(dat_1_5, aes(x=value)) +
  geom_ribbon(aes(ymin=lower, ymax=upper, fill=age), alpha=.15) +
  geom_vline(aes(xintercept=vline), color="black", linetype="dashed", alpha=.4) +
  geom_hline(yintercept=0, color="black", linetype="dashed", alpha=.4) +
  geom_line(aes(y=yhat, color=age), size=.5, alpha=.7) +
  geom_rug(data=filter(dat_1_5, model=="Single model" & sev.only==0 & age=="10 years"), sides="b", alpha=.2, color="red") +
  geom_rug(data=filter(dat_1_5, model=="Single model" & sev.only==1 & age=="10 years"), sides="t", alpha=.2, color="red") +
  geom_rug(data=filter(dat_1_5, model=="Global model" & sev.only==0 & age=="25 years"), sides="b", alpha=.2, color="blue") +
  geom_rug(data=filter(dat_1_5, model=="Global model" & sev.only==1 & age=="25 years"), sides="t", alpha=.2, color="blue") +
  scale_color_manual(values=c("red","blue"), labels=c("Children","Adults")) +
  scale_fill_manual(values=c("red","blue"), labels=c("Children","Adults")) +
  ylab("Odds ratio") +
  theme(legend.position="top", legend.title=element_blank(), axis.title.x=element_blank(), axis.text.x=element_text(angle=45)) +
  facet_grid_sc(cols=vars(biomarker), rows=vars(model), 
                scales=list(x=scales_x, y=scales_y1),
                labeller = as_labeller(c(models, biomarkers)))

# Plot for the last 5 biomarkers
p.6.10 <- ggplot(dat_6_10, aes(x=value)) +
  geom_ribbon(aes(ymin=lower, ymax=upper, fill=age), alpha=.15) +
  geom_vline(aes(xintercept=vline), color="black", linetype="dashed", alpha=.4) +
  geom_hline(yintercept=0, color="black", linetype="dashed", alpha=.4) +
  geom_line(aes(y=yhat, color=age), size=.5, alpha=.7) +
  geom_rug(data=filter(dat_6_10, model=="Single model" & sev.only==0 & age=="10 years"), sides="b", alpha=.2, color="red") +
  geom_rug(data=filter(dat_6_10, model=="Single model" & sev.only==1 & age=="10 years"), sides="t", alpha=.2, color="red") +
  geom_rug(data=filter(dat_6_10, model=="Global model" & sev.only==0 & age=="25 years"), sides="b", alpha=.2, color="blue") +
  geom_rug(data=filter(dat_6_10, model=="Global model" & sev.only==1 & age=="25 years"), sides="t", alpha=.2, color="blue") +
  scale_color_manual(values=c("red","blue"), labels=c("Children","Adults")) +
  scale_fill_manual(values=c("red","blue"), labels=c("Children","Adults")) +
  ylab("Odds ratio") +
  theme(legend.position="none", axis.title.x=element_blank(), axis.text.x=element_text(angle=45)) +
  facet_grid_sc(cols=vars(biomarker), rows=vars(model), 
                scales=list(x=scales_x, y=scales_y2),
                labeller = as_labeller(c(models, biomarkers)))

# Merge plots
p61 <- gridExtra::grid.arrange(p.1.5, p.6.10, nrow=2, heights=c(1,.85))
```

![](Elife-ERA-codes_files/figure-gfm/fig%20A6.1-1.png)<!-- -->

``` r
#ggsave(filename="A6_fig1.pdf", p61, width=8.5, height=6.2)
```

########################################################################################### 

#### Appendix 6—figure 2. Results from models for severe dengue or dengue with warning signs endpoint.

``` r
# Set datadist for 'lrm' function [rms]
dd <- datadist(dat1); options(datadist="dd")

# For children
dd$limits["Adjust to","Age"] <- 10

dat_m1 <- get_pred1(out="sev.or.ws", bio="VCAM", age=10, dat=dat1) %>% rename(value=VCAM) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_m2 <- get_pred1(out="sev.or.ws", bio="SDC", age=10, dat=dat1) %>% rename(value=SDC) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_m3 <- get_pred1(out="sev.or.ws", bio="Ang", age=10, dat=dat1) %>% rename(value=Ang) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_m4 <- get_pred1(out="sev.or.ws", bio="IL8", age=10, dat=dat1) %>% rename(value=IL8) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_m5 <- get_pred1(out="sev.or.ws", bio="IP10", age=10, dat=dat1) %>% rename(value=IP10) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_m6 <- get_pred1(out="sev.or.ws", bio="IL1RA", age=10, dat=dat1) %>% rename(value=IL1RA) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_m7 <- get_pred1(out="sev.or.ws", bio="CD163", age=10, dat=dat1) %>% rename(value=CD163) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_m8 <- get_pred1(out="sev.or.ws", bio="TREM", age=10, dat=dat1) %>% rename(value=TREM) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_m9 <- get_pred1(out="sev.or.ws", bio="Fer", age=10, dat=dat1) %>% rename(value=Fer) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_m10 <- get_pred1(out="sev.or.ws", bio="CRP", age=10, dat=dat1) %>% rename(value=CRP) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))

dat_mg1 <- get_pred2(out="sev.or.ws", bio="VCAM", age=10, dat=dat1) %>% rename(value=VCAM) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_mg2 <- get_pred2(out="sev.or.ws", bio="SDC", age=10, dat=dat1) %>% rename(value=SDC) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_mg3 <- get_pred2(out="sev.or.ws", bio="Ang", age=10, dat=dat1) %>% rename(value=Ang) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_mg4 <- get_pred2(out="sev.or.ws", bio="IL8", age=10, dat=dat1) %>% rename(value=IL8) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_mg5 <- get_pred2(out="sev.or.ws", bio="IP10", age=10, dat=dat1) %>% rename(value=IP10) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_mg6 <- get_pred2(out="sev.or.ws", bio="IL1RA", age=10, dat=dat1) %>% rename(value=IL1RA) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_mg7 <- get_pred2(out="sev.or.ws", bio="CD163", age=10, dat=dat1) %>% rename(value=CD163) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_mg8 <- get_pred2(out="sev.or.ws", bio="TREM", age=10, dat=dat1) %>% rename(value=TREM) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_mg9 <- get_pred2(out="sev.or.ws", bio="Fer", age=10, dat=dat1) %>% rename(value=Fer) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_mg10 <- get_pred2(out="sev.or.ws", bio="CRP", age=10, dat=dat1) %>% rename(value=CRP) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))

dat_p1 <- rbind(dat_m1,dat_mg1, dat_m2,dat_mg2, dat_m3,dat_mg3, dat_m4,dat_mg4, dat_m5,dat_mg5, dat_m6,dat_mg6, dat_m7,dat_mg7, dat_m8,dat_mg8, dat_m9,dat_mg9, dat_m10,dat_mg10) %>%
  mutate(age = "10 years")

# For adults
dd$limits["Adjust to","Age"] <- 25

dat_m1 <- get_pred1(out="sev.or.ws", bio="VCAM", age=25, dat=dat1) %>% rename(value=VCAM) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_m2 <- get_pred1(out="sev.or.ws", bio="SDC", age=25, dat=dat1) %>% rename(value=SDC) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_m3 <- get_pred1(out="sev.or.ws", bio="Ang", age=25, dat=dat1) %>% rename(value=Ang) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_m4 <- get_pred1(out="sev.or.ws", bio="IL8", age=25, dat=dat1) %>% rename(value=IL8) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_m5 <- get_pred1(out="sev.or.ws", bio="IP10", age=25, dat=dat1) %>% rename(value=IP10) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_m6 <- get_pred1(out="sev.or.ws", bio="IL1RA", age=25, dat=dat1) %>% rename(value=IL1RA) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_m7 <- get_pred1(out="sev.or.ws", bio="CD163", age=25, dat=dat1) %>% rename(value=CD163) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_m8 <- get_pred1(out="sev.or.ws", bio="TREM", age=25, dat=dat1) %>% rename(value=TREM) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_m9 <- get_pred1(out="sev.or.ws", bio="Fer", age=25, dat=dat1) %>% rename(value=Fer) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_m10 <- get_pred1(out="sev.or.ws", bio="CRP", age=25, dat=dat1) %>% rename(value=CRP) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))

dat_mg1 <- get_pred2(out="sev.or.ws", bio="VCAM", age=25, dat=dat1) %>% rename(value=VCAM) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_mg2 <- get_pred2(out="sev.or.ws", bio="SDC", age=25, dat=dat1) %>% rename(value=SDC) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_mg3 <- get_pred2(out="sev.or.ws", bio="Ang", age=25, dat=dat1) %>% rename(value=Ang) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_mg4 <- get_pred2(out="sev.or.ws", bio="IL8", age=25, dat=dat1) %>% rename(value=IL8) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_mg5 <- get_pred2(out="sev.or.ws", bio="IP10", age=25, dat=dat1) %>% rename(value=IP10) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_mg6 <- get_pred2(out="sev.or.ws", bio="IL1RA", age=25, dat=dat1) %>% rename(value=IL1RA) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_mg7 <- get_pred2(out="sev.or.ws", bio="CD163", age=25, dat=dat1) %>% rename(value=CD163) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_mg8 <- get_pred2(out="sev.or.ws", bio="TREM", age=25, dat=dat1) %>% rename(value=TREM) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_mg9 <- get_pred2(out="sev.or.ws", bio="Fer", age=25, dat=dat1) %>% rename(value=Fer) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_mg10 <- get_pred2(out="sev.or.ws", bio="CRP", age=25, dat=dat1) %>% rename(value=CRP) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))

dat_p2 <- rbind(dat_m1,dat_mg1, dat_m2,dat_mg2, dat_m3,dat_mg3, dat_m4,dat_mg4, dat_m5,dat_mg5, dat_m6,dat_mg6, dat_m7,dat_mg7, dat_m8,dat_mg8, dat_m9,dat_mg9, dat_m10,dat_mg10) %>%
  mutate(age = "25 years")

# Merge data for children and adults
vis3 <- rbind(dat_p1, dat_p2) %>% mutate(outcome = "outcome 3")

# Merge data for plots
tmp0 <- dat %>% arrange(age15) %>% select(sev.or.ws)
tmp <- data.frame(sev.or.ws = rep(tmp0$sev.or.ws, 20))

tmp_p1 <- vis3 %>%
  mutate(value1 = 2^value,
         age = factor(age, levels=c("10 years", "25 years")),
         gr = ifelse(model=="Single model" & age=="10 years", 1,
                     ifelse(model=="Single model" & age=="25 years", 2,
                            ifelse(model=="Global model" & age=="10 years", 3, 4))),
         gr = factor(gr, levels=c(1:4), labels=c("Single children", "Single adults", "Global children", "Global adults")),
         model = factor(model, levels=c("Single model", "Global model"))) %>%
  bind_cols(., tmp)

vline <- tmp_p1 %>%
  filter(model=="Single model") %>%
  filter(yhat==0) %>%
  select(biomarker, value, value1) %>%
  rename(vline=value, vline1=value1)

dat_p <- left_join(tmp_p1, vline, by="biomarker")

# Alternative data to limit to 5th - 95th of the values
dat_alt <- dat_p %>%
  group_by(biomarker, model) %>%
  mutate(up = quantile(value, .95),
         lo = quantile(value, .05)) %>%
  ungroup() %>%
  mutate(is.outlier = value<lo | value>up | value<0,
         value = ifelse(is.outlier, NA, value),
         value1 = ifelse(is.outlier, NA, value1),
         yhat = ifelse(is.outlier, NA, yhat),
         lower = ifelse(is.outlier, NA, lower),
         upper = ifelse(is.outlier, NA, upper))

dat_1_5 <- dat_alt %>% 
  filter(biomarker %in% c("VCAM", "SDC", "Ang", "IL8", "IP10")) %>%
  mutate(biomarker = factor(biomarker, levels=c("VCAM", "SDC", "Ang", "IL8", "IP10")))

dat_6_10 <- dat_alt %>% 
  filter(biomarker %in% c("IL1RA", "CD163", "TREM", "Fer", "CRP")) %>%
  mutate(biomarker = factor(biomarker, levels=c("IL1RA", "CD163", "TREM", "Fer", "CRP")))

# Modify facets' scales
#library(facetscales)

xVCAM <- c(1,4,15,60,250,1000,4000)
xSDC <- c(1400,2000,2800,4000,5600)
xAng <- c(50,100,250,500,1000,2000)
xIL8 <- c(5,7,10,14,20,28,40)
xIP10 <- c(25,100,400,1600,6400)
xIL1RA <- c(1000,2000,4000,8000,16000)
xCD163 <- c(75,150,300,600)
xTREM <- c(35,50,70,100,140,200)
xFer <- c(50,100,200,400,800)
xCRP <- c(2.5,5,10,20,40,80)

scales_x <- list(
  `VCAM` = scale_x_continuous(breaks = log2(xVCAM), labels = xVCAM),
  `SDC` = scale_x_continuous(breaks = log2(xSDC), labels = xSDC),
  `Ang` = scale_x_continuous(breaks = log2(xAng), labels = xAng), 
  `IL8` = scale_x_continuous(breaks = log2(xIL8), labels = xIL8), 
  `IP10` = scale_x_continuous(breaks = log2(xIP10), labels = xIP10), 
  `IL1RA` = scale_x_continuous(breaks = log2(xIL1RA), labels = xIL1RA), 
  `CD163` = scale_x_continuous(breaks = log2(xCD163), labels = xCD163), 
  `TREM` = scale_x_continuous(breaks = log2(xTREM), labels = xTREM), 
  `Fer` = scale_x_continuous(breaks = log2(xFer), labels = xFer), 
  `CRP` = scale_x_continuous(breaks = log2(xCRP), labels = xCRP)
)

ybreak1 <- c(.25, .5, 1, 2, 4)
ybreak2 <- c(.25, .5, 1, 2, 4)

scales_y1 <- list(
  `Single model` = scale_y_continuous(limits = c(log(.24), log(4.01)), breaks = log(ybreak1), labels = ybreak1),
  `Global model` = scale_y_continuous(limits = c(log(.24), log(4.01)), breaks = log(ybreak1), labels = ybreak1)
)

scales_y2 <- list(
  `Single model` = scale_y_continuous(limits = c(log(.24), log(4.01)), breaks = log(ybreak2), labels = ybreak2),
  `Global model` = scale_y_continuous(limits = c(log(.24), log(4.01)), breaks = log(ybreak2), labels = ybreak2)
)

# Set facets' names
models <- c(`Single model` = "Single model",
            `Global model` = "Global model")

biomarkers <- c(`VCAM` = "VCAM-1 (ng/ml)",
                `SDC`= "SDC-1 (pg/ml)",
                `Ang` = "Ang-2 (pg/ml)",
                `IL8` = "IL-8 (pg/ml)",
                `IP10` = "IP-10 (pg/ml)",
                `IL1RA` = "IL-1RA (pg/ml)",
                `CD163` = "sCD163 (ng/ml)",
                `TREM` = "sTREM-1 (pg/ml)",
                `Fer` = "Ferritin (ng/ml)",
                `CRP` = "CRP (mg/l)")

# Plot for the first 5 biomarkers
p.1.5 <- ggplot(dat_1_5, aes(x=value)) +
  geom_ribbon(aes(ymin=lower, ymax=upper, fill=age), alpha=.15) +
  geom_vline(aes(xintercept=vline), color="black", linetype="dashed", alpha=.4) +
  geom_hline(yintercept=0, color="black", linetype="dashed", alpha=.4) +
  geom_line(aes(y=yhat, color=age), size=.5, alpha=.7) +
  geom_rug(data=filter(dat_1_5, model=="Single model" & sev.or.ws==0 & age=="10 years"), sides="b", alpha=.2, color="red") +
  geom_rug(data=filter(dat_1_5, model=="Single model" & sev.or.ws==1 & age=="10 years"), sides="t", alpha=.2, color="red") +
  geom_rug(data=filter(dat_1_5, model=="Global model" & sev.or.ws==0 & age=="25 years"), sides="b", alpha=.2, color="blue") +
  geom_rug(data=filter(dat_1_5, model=="Global model" & sev.or.ws==1 & age=="25 years"), sides="t", alpha=.2, color="blue") +
  scale_color_manual(values=c("red","blue"), labels=c("Children","Adults")) +
  scale_fill_manual(values=c("red","blue"), labels=c("Children","Adults")) +
  ylab("Odds ratio") +
  theme(legend.position="top", legend.title=element_blank(), axis.title.x=element_blank(), axis.text.x=element_text(angle=45)) +
  facet_grid_sc(cols=vars(biomarker), rows=vars(model), 
                scales=list(x=scales_x, y=scales_y1),
                labeller = as_labeller(c(models, biomarkers)))

# Plot for the last 5 biomarkers
p.6.10 <- ggplot(dat_6_10, aes(x=value)) +
  geom_ribbon(aes(ymin=lower, ymax=upper, fill=age), alpha=.15) +
  geom_vline(aes(xintercept=vline), color="black", linetype="dashed", alpha=.4) +
  geom_hline(yintercept=0, color="black", linetype="dashed", alpha=.4) +
  geom_line(aes(y=yhat, color=age), size=.5, alpha=.7) +
  geom_rug(data=filter(dat_6_10, model=="Single model" & sev.or.ws==0 & age=="10 years"), sides="b", alpha=.2, color="red") +
  geom_rug(data=filter(dat_6_10, model=="Single model" & sev.or.ws==1 & age=="10 years"), sides="t", alpha=.2, color="red") +
  geom_rug(data=filter(dat_6_10, model=="Global model" & sev.or.ws==0 & age=="25 years"), sides="b", alpha=.2, color="blue") +
  geom_rug(data=filter(dat_6_10, model=="Global model" & sev.or.ws==1 & age=="25 years"), sides="t", alpha=.2, color="blue") +
  scale_color_manual(values=c("red","blue"), labels=c("Children","Adults")) +
  scale_fill_manual(values=c("red","blue"), labels=c("Children","Adults")) +
  ylab("Odds ratio") +
  theme(legend.position="none", axis.title.x=element_blank(), axis.text.x=element_text(angle=45)) +
  facet_grid_sc(cols=vars(biomarker), rows=vars(model), 
                scales=list(x=scales_x, y=scales_y2),
                labeller = as_labeller(c(models, biomarkers)))

# Merge plots
p62 <- gridExtra::grid.arrange(p.1.5, p.6.10, nrow=2, heights=c(1,.85))
```

![](Elife-ERA-codes_files/figure-gfm/fig%20A6.2-1.png)<!-- -->

``` r
#ggsave(filename="A6_fig2.pdf", p62, width=8.5, height=6.2)
```

########################################################################################### 

#### Appendix 6—figure 3. Results from models for hospitalization endpoint.

``` r
# Set datadist for 'lrm' function [rms]
dd <- datadist(dat1); options(datadist="dd")

# For children
dd$limits["Adjust to","Age"] <- 10

dat_m1 <- get_pred1(out="hospital", bio="VCAM", age=10, dat=dat1) %>% rename(value=VCAM) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_m2 <- get_pred1(out="hospital", bio="SDC", age=10, dat=dat1) %>% rename(value=SDC) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_m3 <- get_pred1(out="hospital", bio="Ang", age=10, dat=dat1) %>% rename(value=Ang) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_m4 <- get_pred1(out="hospital", bio="IL8", age=10, dat=dat1) %>% rename(value=IL8) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_m5 <- get_pred1(out="hospital", bio="IP10", age=10, dat=dat1) %>% rename(value=IP10) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_m6 <- get_pred1(out="hospital", bio="IL1RA", age=10, dat=dat1) %>% rename(value=IL1RA) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_m7 <- get_pred1(out="hospital", bio="CD163", age=10, dat=dat1) %>% rename(value=CD163) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_m8 <- get_pred1(out="hospital", bio="TREM", age=10, dat=dat1) %>% rename(value=TREM) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_m9 <- get_pred1(out="hospital", bio="Fer", age=10, dat=dat1) %>% rename(value=Fer) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_m10 <- get_pred1(out="hospital", bio="CRP", age=10, dat=dat1) %>% rename(value=CRP) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))

dat_mg1 <- get_pred2(out="hospital", bio="VCAM", age=10, dat=dat1) %>% rename(value=VCAM) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_mg2 <- get_pred2(out="hospital", bio="SDC", age=10, dat=dat1) %>% rename(value=SDC) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_mg3 <- get_pred2(out="hospital", bio="Ang", age=10, dat=dat1) %>% rename(value=Ang) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_mg4 <- get_pred2(out="hospital", bio="IL8", age=10, dat=dat1) %>% rename(value=IL8) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_mg5 <- get_pred2(out="hospital", bio="IP10", age=10, dat=dat1) %>% rename(value=IP10) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_mg6 <- get_pred2(out="hospital", bio="IL1RA", age=10, dat=dat1) %>% rename(value=IL1RA) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_mg7 <- get_pred2(out="hospital", bio="CD163", age=10, dat=dat1) %>% rename(value=CD163) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_mg8 <- get_pred2(out="hospital", bio="TREM", age=10, dat=dat1) %>% rename(value=TREM) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_mg9 <- get_pred2(out="hospital", bio="Fer", age=10, dat=dat1) %>% rename(value=Fer) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_mg10 <- get_pred2(out="hospital", bio="CRP", age=10, dat=dat1) %>% rename(value=CRP) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))

dat_p1 <- rbind(dat_m1,dat_mg1, dat_m2,dat_mg2, dat_m3,dat_mg3, dat_m4,dat_mg4, dat_m5,dat_mg5, dat_m6,dat_mg6, dat_m7,dat_mg7, dat_m8,dat_mg8, dat_m9,dat_mg9, dat_m10,dat_mg10) %>%
  mutate(age = "10 years")

# For adults
dd$limits["Adjust to","Age"] <- 25

dat_m1 <- get_pred1(out="hospital", bio="VCAM", age=25, dat=dat1) %>% rename(value=VCAM) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_m2 <- get_pred1(out="hospital", bio="SDC", age=25, dat=dat1) %>% rename(value=SDC) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_m3 <- get_pred1(out="hospital", bio="Ang", age=25, dat=dat1) %>% rename(value=Ang) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_m4 <- get_pred1(out="hospital", bio="IL8", age=25, dat=dat1) %>% rename(value=IL8) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_m5 <- get_pred1(out="hospital", bio="IP10", age=25, dat=dat1) %>% rename(value=IP10) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_m6 <- get_pred1(out="hospital", bio="IL1RA", age=25, dat=dat1) %>% rename(value=IL1RA) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_m7 <- get_pred1(out="hospital", bio="CD163", age=25, dat=dat1) %>% rename(value=CD163) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_m8 <- get_pred1(out="hospital", bio="TREM", age=25, dat=dat1) %>% rename(value=TREM) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_m9 <- get_pred1(out="hospital", bio="Fer", age=25, dat=dat1) %>% rename(value=Fer) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_m10 <- get_pred1(out="hospital", bio="CRP", age=25, dat=dat1) %>% rename(value=CRP) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))

dat_mg1 <- get_pred2(out="hospital", bio="VCAM", age=25, dat=dat1) %>% rename(value=VCAM) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_mg2 <- get_pred2(out="hospital", bio="SDC", age=25, dat=dat1) %>% rename(value=SDC) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_mg3 <- get_pred2(out="hospital", bio="Ang", age=25, dat=dat1) %>% rename(value=Ang) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_mg4 <- get_pred2(out="hospital", bio="IL8", age=25, dat=dat1) %>% rename(value=IL8) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_mg5 <- get_pred2(out="hospital", bio="IP10", age=25, dat=dat1) %>% rename(value=IP10) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_mg6 <- get_pred2(out="hospital", bio="IL1RA", age=25, dat=dat1) %>% rename(value=IL1RA) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_mg7 <- get_pred2(out="hospital", bio="CD163", age=25, dat=dat1) %>% rename(value=CD163) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_mg8 <- get_pred2(out="hospital", bio="TREM", age=25, dat=dat1) %>% rename(value=TREM) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_mg9 <- get_pred2(out="hospital", bio="Fer", age=25, dat=dat1) %>% rename(value=Fer) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_mg10 <- get_pred2(out="hospital", bio="CRP", age=25, dat=dat1) %>% rename(value=CRP) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))

dat_p2 <- rbind(dat_m1,dat_mg1, dat_m2,dat_mg2, dat_m3,dat_mg3, dat_m4,dat_mg4, dat_m5,dat_mg5, dat_m6,dat_mg6, dat_m7,dat_mg7, dat_m8,dat_mg8, dat_m9,dat_mg9, dat_m10,dat_mg10) %>%
  mutate(age = "25 years")

# Merge data for children and adults
vis4 <- rbind(dat_p1, dat_p2) %>% mutate(outcome = "outcome 4")

# Merge data for plots
tmp0 <- dat %>% arrange(age15) %>% select(hospital)
tmp <- data.frame(hospital = rep(tmp0$hospital, 20))

tmp_p1 <- vis4 %>%
  mutate(value1 = 2^value,
         age = factor(age, levels=c("10 years", "25 years")),
         gr = ifelse(model=="Single model" & age=="10 years", 1,
                     ifelse(model=="Single model" & age=="25 years", 2,
                            ifelse(model=="Global model" & age=="10 years", 3, 4))),
         gr = factor(gr, levels=c(1:4), labels=c("Single children", "Single adults", "Global children", "Global adults")),
         model = factor(model, levels=c("Single model", "Global model"))) %>%
  bind_cols(., tmp)

vline <- tmp_p1 %>%
  filter(model=="Single model") %>%
  filter(yhat==0) %>%
  select(biomarker, value, value1) %>%
  rename(vline=value, vline1=value1)

dat_p <- left_join(tmp_p1, vline, by="biomarker")

# Alternative data to limit to 5th - 95th of the values
dat_alt <- dat_p %>%
  group_by(biomarker, model) %>%
  mutate(up = quantile(value, .95),
         lo = quantile(value, .05)) %>%
  ungroup() %>%
  mutate(is.outlier = value<lo | value>up | value<0,
         value = ifelse(is.outlier, NA, value),
         value1 = ifelse(is.outlier, NA, value1),
         yhat = ifelse(is.outlier, NA, yhat),
         lower = ifelse(is.outlier, NA, lower),
         upper = ifelse(is.outlier, NA, upper))

dat_1_5 <- dat_alt %>% 
  filter(biomarker %in% c("VCAM", "SDC", "Ang", "IL8", "IP10")) %>%
  mutate(biomarker = factor(biomarker, levels=c("VCAM", "SDC", "Ang", "IL8", "IP10")))

dat_6_10 <- dat_alt %>% 
  filter(biomarker %in% c("IL1RA", "CD163", "TREM", "Fer", "CRP")) %>%
  mutate(biomarker = factor(biomarker, levels=c("IL1RA", "CD163", "TREM", "Fer", "CRP")))

# Modify facets' scales
#library(facetscales)

xVCAM <- c(1,4,15,60,250,1000,4000)
xSDC <- c(1400,2000,2800,4000,5600)
xAng <- c(50,100,250,500,1000,2000)
xIL8 <- c(5,7,10,14,20,28,40)
xIP10 <- c(25,100,400,1600,6400)
xIL1RA <- c(1000,2000,4000,8000,16000)
xCD163 <- c(75,150,300,600)
xTREM <- c(35,50,70,100,140,200)
xFer <- c(50,100,200,400,800)
xCRP <- c(2.5,5,10,20,40,80)

scales_x <- list(
  `VCAM` = scale_x_continuous(breaks = log2(xVCAM), labels = xVCAM),
  `SDC` = scale_x_continuous(breaks = log2(xSDC), labels = xSDC),
  `Ang` = scale_x_continuous(breaks = log2(xAng), labels = xAng), 
  `IL8` = scale_x_continuous(breaks = log2(xIL8), labels = xIL8), 
  `IP10` = scale_x_continuous(breaks = log2(xIP10), labels = xIP10), 
  `IL1RA` = scale_x_continuous(breaks = log2(xIL1RA), labels = xIL1RA), 
  `CD163` = scale_x_continuous(breaks = log2(xCD163), labels = xCD163), 
  `TREM` = scale_x_continuous(breaks = log2(xTREM), labels = xTREM), 
  `Fer` = scale_x_continuous(breaks = log2(xFer), labels = xFer), 
  `CRP` = scale_x_continuous(breaks = log2(xCRP), labels = xCRP)
)

ybreak1 <- c(.125, .25, .5, 1, 2, 4, 8)
ybreak2 <- c(.125, .25, .5, 1, 2, 4, 8)

scales_y1 <- list(
  `Single model` = scale_y_continuous(limits = c(log(.124), log(8.01)), breaks = log(ybreak1), labels = ybreak1),
  `Global model` = scale_y_continuous(limits = c(log(.124), log(8.01)), breaks = log(ybreak1), labels = ybreak1)
)

scales_y2 <- list(
  `Single model` = scale_y_continuous(limits = c(log(.124), log(8.01)), breaks = log(ybreak2), labels = ybreak2),
  `Global model` = scale_y_continuous(limits = c(log(.124), log(8.01)), breaks = log(ybreak2), labels = ybreak2)
)

# Set facets' names
models <- c(`Single model` = "Single model",
            `Global model` = "Global model")

biomarkers <- c(`VCAM` = "VCAM-1 (ng/ml)",
                `SDC`= "SDC-1 (pg/ml)",
                `Ang` = "Ang-2 (pg/ml)",
                `IL8` = "IL-8 (pg/ml)",
                `IP10` = "IP-10 (pg/ml)",
                `IL1RA` = "IL-1RA (pg/ml)",
                `CD163` = "sCD163 (ng/ml)",
                `TREM` = "sTREM-1 (pg/ml)",
                `Fer` = "Ferritin (ng/ml)",
                `CRP` = "CRP (mg/l)")

# Plot for the first 5 biomarkers
p.1.5 <- ggplot(dat_1_5, aes(x=value)) +
  geom_ribbon(aes(ymin=lower, ymax=upper, fill=age), alpha=.15) +
  geom_vline(aes(xintercept=vline), color="black", linetype="dashed", alpha=.4) +
  geom_hline(yintercept=0, color="black", linetype="dashed", alpha=.4) +
  geom_line(aes(y=yhat, color=age), size=.5, alpha=.7) +
  geom_rug(data=filter(dat_1_5, model=="Single model" & hospital==0 & age=="10 years"), sides="b", alpha=.2, color="red") +
  geom_rug(data=filter(dat_1_5, model=="Single model" & hospital==1 & age=="10 years"), sides="t", alpha=.2, color="red") +
  geom_rug(data=filter(dat_1_5, model=="Global model" & hospital==0 & age=="25 years"), sides="b", alpha=.2, color="blue") +
  geom_rug(data=filter(dat_1_5, model=="Global model" & hospital==1 & age=="25 years"), sides="t", alpha=.2, color="blue") +
  scale_color_manual(values=c("red","blue"), labels=c("Children","Adults")) +
  scale_fill_manual(values=c("red","blue"), labels=c("Children","Adults")) +
  ylab("Odds ratio") +
  theme(legend.position="top", legend.title=element_blank(), axis.title.x=element_blank(), axis.text.x=element_text(angle=45)) +
  facet_grid_sc(cols=vars(biomarker), rows=vars(model), 
                scales=list(x=scales_x, y=scales_y1),
                labeller = as_labeller(c(models, biomarkers)))

# Plot for the last 5 biomarkers
p.6.10 <- ggplot(dat_6_10, aes(x=value)) +
  geom_ribbon(aes(ymin=lower, ymax=upper, fill=age), alpha=.15) +
  geom_vline(aes(xintercept=vline), color="black", linetype="dashed", alpha=.4) +
  geom_hline(yintercept=0, color="black", linetype="dashed", alpha=.4) +
  geom_line(aes(y=yhat, color=age), size=.5, alpha=.7) +
  geom_rug(data=filter(dat_6_10, model=="Single model" & hospital==0 & age=="10 years"), sides="b", alpha=.2, color="red") +
  geom_rug(data=filter(dat_6_10, model=="Single model" & hospital==1 & age=="10 years"), sides="t", alpha=.2, color="red") +
  geom_rug(data=filter(dat_6_10, model=="Global model" & hospital==0 & age=="25 years"), sides="b", alpha=.2, color="blue") +
  geom_rug(data=filter(dat_6_10, model=="Global model" & hospital==1 & age=="25 years"), sides="t", alpha=.2, color="blue") +
  scale_color_manual(values=c("red","blue"), labels=c("Children","Adults")) +
  scale_fill_manual(values=c("red","blue"), labels=c("Children","Adults")) +
  ylab("Odds ratio") +
  theme(legend.position="none", axis.title.x=element_blank(), axis.text.x=element_text(angle=45)) +
  facet_grid_sc(cols=vars(biomarker), rows=vars(model), 
                scales=list(x=scales_x, y=scales_y2),
                labeller = as_labeller(c(models, biomarkers)))

# Merge plots
p63 <- gridExtra::grid.arrange(p.1.5, p.6.10, nrow=2, heights=c(1,.85))
```

![](Elife-ERA-codes_files/figure-gfm/fig%20A6.3-1.png)<!-- -->

``` r
#ggsave(filename="A6_fig3.pdf", p63, width=8.5, height=6.2)
```

########################################################################################### 

#### Appendix 6—table 1. Results from models for severe dengue endpoint.

``` r
# Set datadist for 'lrm' function [rms]
dd <- datadist(dat1); options(datadist="dd")

# Get ORs and 95% CIs
tmp1 <- data.frame(
  sort1 = rep(c(1:10), 2),
  sort2 = c(rep(2,10), rep(3,10)),
  bio = rep(c("VCAM", "SDC", "Ang", "IL8", "IP10", "IL1RA", "CD163", "TREM", "Fer", "CRP"), 2),
  ref1 = c(ref0-1, ref0),
  ref2 = c(ref0, ref0+1)
) %>%
  arrange(sort1, sort2) %>%
  group_by(sort1, sort2) %>%
  do(cbind(.,
           # Children - single models
           or1c = get_est1(out="sev.only", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="OR", age=10, dat=dat1),
           lo1c = get_est1(out="sev.only", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="loCI", age=10, dat=dat1),
           up1c = get_est1(out="sev.only", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="upCI", age=10, dat=dat1),
           # Children - global model
           or2c = get_est2(out="sev.only", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="OR", age=10, dat=dat1),
           lo2c = get_est2(out="sev.only", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="loCI", age=10, dat=dat1),
           up2c = get_est2(out="sev.only", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="upCI", age=10, dat=dat1),
           # Adults - single models
           or1a = get_est1(out="sev.only", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="OR", age=25, dat=dat1),
           lo1a = get_est1(out="sev.only", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="loCI", age=25, dat=dat1),
           up1a = get_est1(out="sev.only", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="upCI", age=25, dat=dat1),
           # Adults - global model
           or2a = get_est2(out="sev.only", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="OR", age=25, dat=dat1),
           lo2a = get_est2(out="sev.only", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="loCI", age=25, dat=dat1),
           up2a = get_est2(out="sev.only", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="upCI", age=25, dat=dat1))) %>%
  ungroup()
for (i in 6:17) {tmp1[[i]] <- sprintf("%.2f", round(tmp1[[i]],2))}

# Get p-values
tmp2 <- data.frame(
  sort1 = c(1:10),
  sort2 = rep(1,10),
  bio = c("VCAM", "SDC", "Ang", "IL8", "IP10", "IL1RA", "CD163", "TREM", "Fer", "CRP")
) %>%
  group_by(sort1) %>%
  do(cbind(.,
           p.s2 = get_est1(out="sev.only", bio=.$bio, est="p", dat=dat1),
           p.s2_int = get_est1(out="sev.only", bio=.$bio, est="p int", dat=dat1),
           p.g2 = get_est2(out="sev.only", bio=.$bio, est="p", dat=dat1),
           p.g2_int = get_est2(out="sev.only", bio=.$bio, est="p int", dat=dat1))) %>%
  ungroup()
for (i in 4:7) {tmp2[[i]] <- ifelse(tmp2[[i]]<0.001, "<0.001", sprintf("%.3f", round(tmp2[[i]],3)))}  

# Combine results
res2 <- tmp1 %>%
  group_by(bio) %>%
  slice(1) %>%
  left_join(., select(tmp2, -c(sort1, sort2)), by="bio") %>%
  arrange(sort1) %>%
  mutate(bio = as.character(bio),
         or.sc2 = ifelse(is.na(lo1c), NA, paste(or1c, " (", lo1c, "-", up1c, ")", sep="")), # s: single model; c: children
         or.gc2 = ifelse(is.na(lo2c), NA, paste(or2c, " (", lo2c, "-", up2c, ")", sep="")), # g: global model
         or.sa2 = ifelse(is.na(lo1a), NA, paste(or1a, " (", lo1a, "-", up1a, ")", sep="")), # a: adults
         or.ga2 = ifelse(is.na(lo2a), NA, paste(or2a, " (", lo2a, "-", up2a, ")", sep=""))) %>%
  select(bio, or.sc2, or.sa2, p.s2, p.s2_int, or.gc2, or.ga2, p.g2, p.g2_int)

names(res2) <- c("", "OR (single - children)", "OR (single - adults)", "P overall (single)", "P interaction (single)",
                 "OR (global - children)", "OR (global - adults)", "P overall (global)", "P interaction (global)")

# Report results
knitr::kable(res2)
```

|       | OR (single - children) | OR (single - adults) | P overall (single) | P interaction (single) | OR (global - children) | OR (global - adults) | P overall (global) | P interaction (global) |
|:------|:-----------------------|:---------------------|:-------------------|:-----------------------|:-----------------------|:---------------------|:-------------------|:-----------------------|
| VCAM  | 1.28 (0.77-2.12)       | 2.13 (0.88-5.15)     | 0.236              | 0.202                  | 1.20 (0.35-4.12)       | 3.13 (0.19-50.77)    | 0.723              | 0.491                  |
| SDC   | 1.55 (0.41-5.96)       | 3.26 (0.69-15.52)    | 0.307              | 0.438                  | 1.16 (0.12-11.07)      | 0.87 (0.03-27.93)    | 0.987              | 0.884                  |
| Ang   | 1.42 (0.70-2.88)       | 1.79 (0.97-3.31)     | 0.148              | 0.580                  | 1.20 (0.33-4.31)       | 1.02 (0.33-3.19)     | 0.961              | 0.830                  |
| IL8   | 1.04 (0.47-2.33)       | 2.15 (0.67-6.88)     | 0.420              | 0.338                  | 0.81 (0.28-2.31)       | 2.51 (0.16-38.43)    | 0.741              | 0.446                  |
| IP10  | 1.32 (0.75-2.33)       | 1.64 (0.58-4.62)     | 0.443              | 0.708                  | 0.76 (0.17-3.30)       | 0.22 (0.01-7.12)     | 0.693              | 0.473                  |
| IL1RA | 1.84 (0.81-4.16)       | 1.78 (0.80-3.94)     | 0.098              | 0.956                  | 2.22 (0.51-9.64)       | 3.35 (0.45-24.87)    | 0.386              | 0.696                  |
| CD163 | 1.43 (0.52-3.94)       | 2.43 (0.27-21.70)    | 0.617              | 0.645                  | 0.98 (0.27-3.59)       | 2.06 (0.15-27.47)    | 0.855              | 0.590                  |
| TREM  | 1.20 (0.44-3.28)       | 1.29 (0.55-3.05)     | 0.793              | 0.914                  | 0.93 (0.30-2.88)       | 0.38 (0.02-8.76)     | 0.827              | 0.541                  |
| Fer   | 1.13 (0.59-2.16)       | 1.58 (0.59-4.23)     | 0.659              | 0.501                  | 1.12 (0.55-2.29)       | 1.34 (0.25-7.02)     | 0.926              | 0.819                  |
| CRP   | 1.16 (0.71-1.88)       | 1.77 (0.56-5.62)     | 0.620              | 0.418                  | 0.97 (0.58-1.63)       | 1.55 (0.33-7.32)     | 0.773              | 0.500                  |

########################################################################################### 

#### Appendix 6—table 2. Results from models for severe dengue or dengue with warning signs endpoint.

``` r
# Set datadist for 'lrm' function [rms]
dd <- datadist(dat1); options(datadist="dd")

# Get ORs and 95% CIs
tmp1 <- data.frame(
  sort1 = rep(c(1:10), 2),
  sort2 = c(rep(2,10), rep(3,10)),
  bio = rep(c("VCAM", "SDC", "Ang", "IL8", "IP10", "IL1RA", "CD163", "TREM", "Fer", "CRP"), 2),
  ref1 = c(ref0-1, ref0),
  ref2 = c(ref0, ref0+1)
) %>%
  arrange(sort1, sort2) %>%
  group_by(sort1, sort2) %>%
  do(cbind(.,
           # Children - single models
           or1c = get_est1(out="sev.or.ws", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="OR", age=10, dat=dat1),
           lo1c = get_est1(out="sev.or.ws", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="loCI", age=10, dat=dat1),
           up1c = get_est1(out="sev.or.ws", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="upCI", age=10, dat=dat1),
           # Children - global model
           or2c = get_est2(out="sev.or.ws", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="OR", age=10, dat=dat1),
           lo2c = get_est2(out="sev.or.ws", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="loCI", age=10, dat=dat1),
           up2c = get_est2(out="sev.or.ws", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="upCI", age=10, dat=dat1),
           # Adults - single models
           or1a = get_est1(out="sev.or.ws", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="OR", age=25, dat=dat1),
           lo1a = get_est1(out="sev.or.ws", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="loCI", age=25, dat=dat1),
           up1a = get_est1(out="sev.or.ws", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="upCI", age=25, dat=dat1),
           # Adults - global model
           or2a = get_est2(out="sev.or.ws", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="OR", age=25, dat=dat1),
           lo2a = get_est2(out="sev.or.ws", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="loCI", age=25, dat=dat1),
           up2a = get_est2(out="sev.or.ws", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="upCI", age=25, dat=dat1))) %>%
  ungroup()
for (i in 6:17) {tmp1[[i]] <- sprintf("%.2f", round(tmp1[[i]],2))}

# Get p-values
tmp2 <- data.frame(
  sort1 = c(1:10),
  sort2 = rep(1,10),
  bio = c("VCAM", "SDC", "Ang", "IL8", "IP10", "IL1RA", "CD163", "TREM", "Fer", "CRP")
) %>%
  group_by(sort1) %>%
  do(cbind(.,
           p.s3 = get_est1(out="sev.or.ws", bio=.$bio, est="p", dat=dat1),
           p.s3_int = get_est1(out="sev.or.ws", bio=.$bio, est="p int", dat=dat1),
           p.g3 = get_est2(out="sev.or.ws", bio=.$bio, est="p", dat=dat1),
           p.g3_int = get_est2(out="sev.or.ws", bio=.$bio, est="p int", dat=dat1))) %>%
  ungroup()
for (i in 4:7) {tmp2[[i]] <- ifelse(tmp2[[i]]<0.001, "<0.001", sprintf("%.3f", round(tmp2[[i]],3)))}  

# Combine results
res3 <- bind_rows(tmp1, tmp2) %>%
  arrange(sort1, sort2) %>%
  mutate(bio = ifelse(!is.na(ref1), paste(" -", round(2^ref2,0), "vs", round(2^ref1,0), sep=" "), as.character(bio)),
         or.sc3 = ifelse(is.na(lo1c), NA, paste(or1c, " (", lo1c, "-", up1c, ")", sep="")), # s: single model; c: children
         or.gc3 = ifelse(is.na(lo2c), NA, paste(or2c, " (", lo2c, "-", up2c, ")", sep="")), # g: global model
         or.sa3 = ifelse(is.na(lo1a), NA, paste(or1a, " (", lo1a, "-", up1a, ")", sep="")), # a: adults
         or.ga3 = ifelse(is.na(lo2a), NA, paste(or2a, " (", lo2a, "-", up2a, ")", sep=""))) %>%
  select(bio, or.sc3, or.sa3, p.s3, p.s3_int, or.gc3, or.ga3, p.g3, p.g3_int)

names(res3) <- c("", "OR (single - children)", "OR (single - adults)", "P overall (single)", "P interaction (single)",
                 "OR (global - children)", "OR (global - adults)", "P overall (global)", "P interaction (global)")

# Report results
knitr::kable(res3)
```

|                  | OR (single - children) | OR (single - adults) | P overall (single) | P interaction (single) | OR (global - children) | OR (global - adults) | P overall (global) | P interaction (global) |
|:-----------------|:-----------------------|:---------------------|:-------------------|:-----------------------|:-----------------------|:---------------------|:-------------------|:-----------------------|
| VCAM             |                        |                      | 0.025              | 0.374                  |                        |                      | 0.469              | 0.763                  |
| \- 1636 vs 818   | 1.06 (0.93-1.22)       | 1.11 (0.92-1.33)     |                    |                        | 0.82 (0.66-1.02)       | 0.92 (0.66-1.30)     |                    |                        |
| \- 3272 vs 1636  | 1.06 (0.87-1.29)       | 1.13 (0.88-1.46)     |                    |                        | 0.74 (0.55-1.01)       | 0.88 (0.56-1.40)     |                    |                        |
| SDC              |                        |                      | 0.032              | 0.363                  |                        |                      | 0.116              | 0.773                  |
| \- 2519 vs 1260  | 1.34 (0.86-2.09)       | 1.15 (0.56-2.33)     |                    |                        | 1.57 (0.89-2.75)       | 1.89 (0.71-5.03)     |                    |                        |
| \- 5039 vs 2519  | 1.40 (0.90-2.17)       | 2.00 (0.79-5.08)     |                    |                        | 1.78 (1.00-3.18)       | 2.88 (0.71-11.69)    |                    |                        |
| Ang              |                        |                      | 0.008              | 0.637                  |                        |                      | 0.009              | 0.011                  |
| \- 1204 vs 602   | 1.23 (1.05-1.44)       | 1.12 (0.93-1.36)     |                    |                        | 1.37 (1.07-1.75)       | 0.95 (0.69-1.30)     |                    |                        |
| \- 2409 vs 1204  | 1.34 (0.95-1.88)       | 1.22 (0.82-1.80)     |                    |                        | 1.42 (0.91-2.21)       | 0.85 (0.48-1.51)     |                    |                        |
| IL8              |                        |                      | 0.040              | 0.020                  |                        |                      | 0.053              | 0.030                  |
| \- 14 vs 7       | 1.05 (0.87-1.27)       | 0.97 (0.67-1.41)     |                    |                        | 0.96 (0.77-1.21)       | 0.94 (0.62-1.42)     |                    |                        |
| \- 28 vs 14      | 1.01 (0.85-1.19)       | 1.45 (0.94-2.22)     |                    |                        | 0.78 (0.62-0.99)       | 1.48 (0.89-2.44)     |                    |                        |
| IP10             |                        |                      | &lt;0.001          | 0.176                  |                        |                      | 0.059              | 0.537                  |
| \- 3093 vs 1546  | 1.26 (1.09-1.44)       | 1.44 (1.16-1.80)     |                    |                        | 1.16 (0.87-1.55)       | 1.30 (0.81-2.09)     |                    |                        |
| \- 6186 vs 3093  | 1.39 (1.13-1.71)       | 1.75 (1.26-2.44)     |                    |                        | 1.33 (0.89-1.99)       | 1.49 (0.78-2.87)     |                    |                        |
| IL1RA            |                        |                      | 0.005              | 0.381                  |                        |                      | 0.425              | 0.955                  |
| \- 6434 vs 3217  | 1.17 (1.06-1.30)       | 1.15 (0.97-1.36)     |                    |                        | 1.24 (1.01-1.52)       | 1.14 (0.80-1.63)     |                    |                        |
| \- 12868 vs 6434 | 1.37 (1.06-1.77)       | 1.70 (1.13-2.55)     |                    |                        | 1.38 (0.93-2.05)       | 1.45 (0.80-2.64)     |                    |                        |
| CD163            |                        |                      | 0.854              | 0.719                  |                        |                      | 0.193              | 0.419                  |
| \- 295 vs 147    | 0.96 (0.84-1.08)       | 1.03 (0.85-1.25)     |                    |                        | 0.85 (0.71-1.03)       | 1.12 (0.81-1.55)     |                    |                        |
| \- 589 vs 295    | 0.91 (0.69-1.22)       | 1.03 (0.61-1.73)     |                    |                        | 0.71 (0.50-1.00)       | 0.80 (0.44-1.47)     |                    |                        |
| TREM             |                        |                      | 0.221              | 0.472                  |                        |                      | 0.002              | 0.132                  |
| \- 85 vs 42      | 0.75 (0.56-1.00)       | 0.69 (0.41-1.16)     |                    |                        | 0.48 (0.32-0.73)       | 0.64 (0.36-1.16)     |                    |                        |
| \- 169 vs 85     | 1.04 (0.84-1.28)       | 0.78 (0.57-1.07)     |                    |                        | 0.87 (0.69-1.10)       | 0.62 (0.39-1.00)     |                    |                        |
| Fer              |                        |                      | 0.034              | 0.258                  |                        |                      | 0.024              | 0.075                  |
| \- 243 vs 122    | 1.13 (1.02-1.26)       | 1.01 (0.85-1.19)     |                    |                        | 1.17 (1.01-1.35)       | 0.96 (0.76-1.23)     |                    |                        |
| \- 487 vs 243    | 0.96 (0.77-1.20)       | 0.76 (0.53-1.08)     |                    |                        | 0.94 (0.71-1.25)       | 0.67 (0.42-1.07)     |                    |                        |
| CRP              |                        |                      | 0.747              | 0.622                  |                        |                      | 0.662              | 0.448                  |
| \- 28 vs 14      | 1.03 (0.96-1.12)       | 0.97 (0.79-1.20)     |                    |                        | 0.99 (0.89-1.10)       | 0.84 (0.59-1.20)     |                    |                        |
| \- 56 vs 28      | 1.02 (0.87-1.19)       | 1.19 (0.92-1.54)     |                    |                        | 0.96 (0.79-1.17)       | 1.06 (0.75-1.49)     |                    |                        |

########################################################################################### 

#### Appendix 6—table 3. Results from models for hospitalization endpoint.

``` r
# Set datadist for 'lrm' function [rms]
dd <- datadist(dat1); options(datadist="dd")

# Get ORs and 95% CIs
tmp1 <- data.frame(
  sort1 = rep(c(1:10), 2),
  sort2 = c(rep(2,10), rep(3,10)),
  bio = rep(c("VCAM", "SDC", "Ang", "IL8", "IP10", "IL1RA", "CD163", "TREM", "Fer", "CRP"), 2),
  ref1 = c(ref0-1, ref0),
  ref2 = c(ref0, ref0+1)
) %>%
  arrange(sort1, sort2) %>%
  group_by(sort1, sort2) %>%
  do(cbind(.,
           # Children - single models
           or1c = get_est1(out="hospital", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="OR", age=10, dat=dat1),
           lo1c = get_est1(out="hospital", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="loCI", age=10, dat=dat1),
           up1c = get_est1(out="hospital", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="upCI", age=10, dat=dat1),
           # Children - global model
           or2c = get_est2(out="hospital", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="OR", age=10, dat=dat1),
           lo2c = get_est2(out="hospital", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="loCI", age=10, dat=dat1),
           up2c = get_est2(out="hospital", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="upCI", age=10, dat=dat1),
           # Adults - single models
           or1a = get_est1(out="hospital", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="OR", age=25, dat=dat1),
           lo1a = get_est1(out="hospital", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="loCI", age=25, dat=dat1),
           up1a = get_est1(out="hospital", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="upCI", age=25, dat=dat1),
           # Adults - global model
           or2a = get_est2(out="hospital", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="OR", age=25, dat=dat1),
           lo2a = get_est2(out="hospital", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="loCI", age=25, dat=dat1),
           up2a = get_est2(out="hospital", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="upCI", age=25, dat=dat1))) %>%
  ungroup()
for (i in 6:17) {tmp1[[i]] <- sprintf("%.2f", round(tmp1[[i]],2))}

# Get p-values
tmp2 <- data.frame(
  sort1 = c(1:10),
  sort2 = rep(1,10),
  bio = c("VCAM", "SDC", "Ang", "IL8", "IP10", "IL1RA", "CD163", "TREM", "Fer", "CRP")
) %>%
  group_by(sort1) %>%
  do(cbind(.,
           p.s4 = get_est1(out="hospital", bio=.$bio, est="p", dat=dat1),
           p.s4_int = get_est1(out="hospital", bio=.$bio, est="p int", dat=dat1),
           p.g4 = get_est2(out="hospital", bio=.$bio, est="p", dat=dat1),
           p.g4_int = get_est2(out="hospital", bio=.$bio, est="p int", dat=dat1))) %>%
  ungroup()
for (i in 4:7) {tmp2[[i]] <- ifelse(tmp2[[i]]<0.001, "<0.001", sprintf("%.3f", round(tmp2[[i]],3)))}  

# Combine results
res4 <- bind_rows(tmp1, tmp2) %>%
  arrange(sort1, sort2) %>%
  mutate(bio = ifelse(!is.na(ref1), paste(" -", round(2^ref2,0), "vs", round(2^ref1,0), sep=" "), as.character(bio)),
         or.sc4 = ifelse(is.na(lo1c), NA, paste(or1c, " (", lo1c, "-", up1c, ")", sep="")), # s: single model; c: children
         or.gc4 = ifelse(is.na(lo2c), NA, paste(or2c, " (", lo2c, "-", up2c, ")", sep="")), # g: global model
         or.sa4 = ifelse(is.na(lo1a), NA, paste(or1a, " (", lo1a, "-", up1a, ")", sep="")), # a: adults
         or.ga4 = ifelse(is.na(lo2a), NA, paste(or2a, " (", lo2a, "-", up2a, ")", sep=""))) %>%
  select(bio, or.sc4, or.sa4, p.s4, p.s4_int, or.gc4, or.ga4, p.g4, p.g4_int)

names(res4) <- c("", "OR (single - children)", "OR (single - adults)", "P overall (single)", "P interaction (single)",
                 "OR (global - children)", "OR (global - adults)", "P overall (global)", "P interaction (global)")

# Report results
knitr::kable(res4)
```

|                  | OR (single - children) | OR (single - adults) | P overall (single) | P interaction (single) | OR (global - children) | OR (global - adults) | P overall (global) | P interaction (global) |
|:-----------------|:-----------------------|:---------------------|:-------------------|:-----------------------|:-----------------------|:---------------------|:-------------------|:-----------------------|
| VCAM             |                        |                      | &lt;0.001          | 0.009                  |                        |                      | 0.092              | 0.137                  |
| \- 1636 vs 818   | 1.28 (1.11-1.49)       | 1.28 (1.03-1.60)     |                    |                        | 1.16 (0.90-1.48)       | 1.28 (0.87-1.90)     |                    |                        |
| \- 3272 vs 1636  | 1.42 (1.14-1.76)       | 1.46 (1.08-1.97)     |                    |                        | 1.16 (0.81-1.64)       | 1.41 (0.84-2.37)     |                    |                        |
| SDC              |                        |                      | &lt;0.001          | 0.187                  |                        |                      | 0.006              | 0.406                  |
| \- 2519 vs 1260  | 1.82 (1.01-3.28)       | 0.79 (0.31-2.05)     |                    |                        | 1.70 (0.84-3.41)       | 0.62 (0.16-2.49)     |                    |                        |
| \- 5039 vs 2519  | 2.22 (1.36-3.63)       | 1.81 (0.65-5.07)     |                    |                        | 3.70 (1.90-7.22)       | 1.66 (0.39-7.07)     |                    |                        |
| Ang              |                        |                      | 0.012              | 0.337                  |                        |                      | 0.497              | 0.789                  |
| \- 1204 vs 602   | 1.27 (1.07-1.52)       | 1.21 (0.91-1.62)     |                    |                        | 1.27 (0.92-1.77)       | 1.16 (0.76-1.77)     |                    |                        |
| \- 2409 vs 1204  | 1.58 (1.08-2.32)       | 1.61 (0.86-3.04)     |                    |                        | 1.63 (0.90-2.93)       | 1.45 (0.70-3.01)     |                    |                        |
| IL8              |                        |                      | 0.007              | 0.002                  |                        |                      | 0.024              | 0.021                  |
| \- 14 vs 7       | 1.14 (0.90-1.44)       | 0.73 (0.47-1.15)     |                    |                        | 0.94 (0.66-1.33)       | 0.63 (0.36-1.13)     |                    |                        |
| \- 28 vs 14      | 0.98 (0.81-1.18)       | 1.32 (0.85-2.05)     |                    |                        | 0.70 (0.50-0.97)       | 1.59 (0.81-3.12)     |                    |                        |
| IP10             |                        |                      | 0.002              | 0.242                  |                        |                      | 0.005              | 0.212                  |
| \- 3093 vs 1546  | 1.32 (1.13-1.54)       | 1.10 (0.86-1.40)     |                    |                        | 0.80 (0.56-1.14)       | 0.77 (0.45-1.30)     |                    |                        |
| \- 6186 vs 3093  | 1.50 (1.18-1.90)       | 1.17 (0.81-1.70)     |                    |                        | 0.84 (0.51-1.40)       | 0.66 (0.32-1.35)     |                    |                        |
| IL1RA            |                        |                      | &lt;0.001          | 0.685                  |                        |                      | &lt;0.001          | 0.389                  |
| \- 6434 vs 3217  | 1.31 (1.17-1.46)       | 1.13 (0.94-1.36)     |                    |                        | 2.05 (1.57-2.66)       | 1.18 (0.73-1.89)     |                    |                        |
| \- 12868 vs 6434 | 1.41 (1.11-1.80)       | 1.17 (0.79-1.72)     |                    |                        | 2.25 (1.39-3.64)       | 1.19 (0.62-2.28)     |                    |                        |
| CD163            |                        |                      | 0.208              | 0.722                  |                        |                      | 0.007              | 0.419                  |
| \- 295 vs 147    | 0.90 (0.79-1.04)       | 0.98 (0.76-1.27)     |                    |                        | 0.69 (0.53-0.89)       | 0.96 (0.62-1.50)     |                    |                        |
| \- 589 vs 295    | 0.69 (0.47-1.00)       | 0.91 (0.46-1.82)     |                    |                        | 0.49 (0.30-0.80)       | 0.83 (0.36-1.91)     |                    |                        |
| TREM             |                        |                      | 0.635              | 0.371                  |                        |                      | 0.011              | 0.053                  |
| \- 85 vs 42      | 0.93 (0.67-1.29)       | 1.35 (0.71-2.58)     |                    |                        | 0.53 (0.34-0.81)       | 1.74 (0.67-4.52)     |                    |                        |
| \- 169 vs 85     | 0.89 (0.71-1.13)       | 0.83 (0.54-1.28)     |                    |                        | 0.55 (0.36-0.83)       | 0.57 (0.26-1.29)     |                    |                        |
| Fer              |                        |                      | &lt;0.001          | 0.011                  |                        |                      | 0.129              | 0.117                  |
| \- 243 vs 122    | 1.22 (1.08-1.38)       | 0.92 (0.72-1.18)     |                    |                        | 1.23 (1.02-1.49)       | 0.91 (0.68-1.21)     |                    |                        |
| \- 487 vs 243    | 1.40 (1.10-1.79)       | 0.99 (0.68-1.46)     |                    |                        | 1.25 (0.89-1.75)       | 0.93 (0.54-1.59)     |                    |                        |
| CRP              |                        |                      | 0.379              | 0.190                  |                        |                      | 0.139              | 0.053                  |
| \- 28 vs 14      | 0.97 (0.89-1.06)       | 1.01 (0.83-1.24)     |                    |                        | 0.97 (0.86-1.11)       | 0.95 (0.71-1.27)     |                    |                        |
| \- 56 vs 28      | 0.95 (0.78-1.15)       | 1.35 (1.01-1.81)     |                    |                        | 0.89 (0.69-1.15)       | 1.37 (0.88-2.13)     |                    |                        |

########################################################################################### 

#### Appendix 7—table 1. Model selection frequencies for children.

``` r
# EPV --------------------------------------------------------------
pred <- c("VCAM", "SDC", "Ang", "IL8", "IP10", "IL1RA", "CD163", "TREM", "Fer", "CRP")

# Bootstrap results ------------------------------------------------
## The following codes were modified from Heinze G. et al (https://doi.org/10.1002/bimj.201700067)
bootnum <- 1000
boot_var <- matrix(0, ncol = length(pred) + 1, nrow = bootnum, dimnames = list(NULL, c(pred, "aic")))
load("5a.RData") # Load bootstrap results

# Model frequency
group_cols <- c(colnames(data.frame(boot_var))[1:length(pred)])
boot_modfreq <- data.frame(boot_var) %>%
  mutate(time = 1) %>%
  group_by_(.dots = group_cols) %>%
  summarise(aic_median = round(median(aic),1),
            aic_025 = round(quantile(aic, 0.025),1),
            aic_975 = round(quantile(aic, 0.975),1),
            count = sum(time)) %>%
  ungroup() %>%
  mutate(percent = count/bootnum * 100) %>%
  arrange(desc(count)) %>%
  slice(1:20) %>%
  as.data.frame(.)

out <- boot_modfreq
for (i in 1:(ncol(out)-5)) {
  out[,i] <- ifelse(out[,i]==0, NA, "+")
}

knitr::kable(out)
```

| VCAM | SDC | Ang | IL8 | IP10 | IL1RA | CD163 | TREM | Fer | CRP | aic\_median | aic\_025 | aic\_975 | count | percent |
|:-----|:----|:----|:----|:-----|:------|:------|:-----|:----|:----|------------:|---------:|---------:|------:|--------:|
|      | \+  | \+  | \+  | \+   | \+    |       |      | \+  |     |       457.9 |    413.0 |    500.7 |   134 |    13.4 |
|      | \+  | \+  | \+  | \+   | \+    | \+    |      | \+  |     |       453.5 |    415.8 |    495.6 |   100 |    10.0 |
|      |     | \+  | \+  | \+   | \+    | \+    |      | \+  |     |       459.7 |    401.8 |    500.6 |    55 |     5.5 |
|      |     | \+  | \+  | \+   | \+    |       |      | \+  |     |       459.2 |    420.8 |    490.4 |    54 |     5.4 |
| \+   | \+  | \+  | \+  |      | \+    | \+    |      | \+  |     |       462.1 |    420.2 |    507.5 |    48 |     4.8 |
| \+   | \+  | \+  | \+  | \+   | \+    |       |      | \+  |     |       450.2 |    422.3 |    492.1 |    47 |     4.7 |
| \+   | \+  | \+  | \+  | \+   | \+    | \+    |      | \+  |     |       460.6 |    418.4 |    502.0 |    46 |     4.6 |
| \+   | \+  | \+  | \+  |      | \+    |       |      | \+  |     |       457.6 |    417.7 |    497.6 |    40 |     4.0 |
|      | \+  | \+  | \+  | \+   | \+    |       |      | \+  | \+  |       451.9 |    413.7 |    493.5 |    39 |     3.9 |
|      | \+  | \+  | \+  | \+   | \+    | \+    | \+   | \+  |     |       453.2 |    394.6 |    491.9 |    36 |     3.6 |
|      | \+  | \+  | \+  | \+   | \+    | \+    |      | \+  | \+  |       452.8 |    418.2 |    484.7 |    28 |     2.8 |
|      | \+  | \+  | \+  | \+   | \+    |       | \+   | \+  |     |       441.9 |    407.8 |    486.8 |    23 |     2.3 |
| \+   |     | \+  | \+  | \+   | \+    |       |      | \+  |     |       461.4 |    434.5 |    490.8 |    23 |     2.3 |
| \+   |     | \+  | \+  |      | \+    |       |      | \+  |     |       466.8 |    427.4 |    481.0 |    17 |     1.7 |
|      | \+  | \+  | \+  | \+   | \+    |       | \+   | \+  | \+  |       449.9 |    431.1 |    482.4 |    15 |     1.5 |
| \+   | \+  | \+  | \+  |      | \+    |       |      | \+  | \+  |       449.3 |    423.0 |    486.2 |    14 |     1.4 |
|      |     | \+  | \+  | \+   | \+    |       |      | \+  | \+  |       449.1 |    420.8 |    477.8 |    13 |     1.3 |
| \+   |     | \+  | \+  | \+   | \+    | \+    |      | \+  |     |       474.9 |    442.2 |    488.7 |    12 |     1.2 |
| \+   | \+  | \+  | \+  | \+   | \+    |       | \+   | \+  |     |       451.0 |    416.7 |    471.9 |    12 |     1.2 |
|      |     | \+  | \+  | \+   | \+    | \+    |      | \+  | \+  |       453.9 |    435.2 |    492.2 |    11 |     1.1 |

########################################################################################### 

#### Appendix 7—table 2. Model stability for children.

``` r
# EPV --------------------------------------------------------------
pred <- c("VCAM", "SDC", "Ang", "IL8", "IP10", "IL1RA", "CD163", "TREM", "Fer", "CRP")

# Estimate full model ----------------------------------------------
full_mod <- glm(sev.or.inte ~ VCAM + SDC + Ang + IL8 + IP10 + IL1RA + CD163 + TREM + Fer + CRP, 
                family=binomial, data=dat1c, x=T, y=T)
full_est <- coef(full_mod)
full_se <- coef(summary(full_mod))[, "Std. Error"]
pred_name <- names(full_est)[-1]

# Selected model (for bootstrap results) ---------------------------
sel_m <- dredge(full_mod, rank="AIC")
sel_var <- matrix(0, ncol = length(pred), nrow = 1, dimnames = list(NULL, pred))
sel_var_tmp <- attr(get.models(sel_m, 1)[[1]]$terms, "term.labels")

for (j in 1:ncol(sel_var)) {
  sel_var[1,j] <- ifelse(names(sel_var[1,j]) %in% sel_var_tmp, 1, 0)
}

formula <- paste("sev.or.inte~", paste(names(sel_var[1,][sel_var[1,]==1]), collapse = "+"))
sel_mod <- glm(formula, data = dat1, family = binomial, x = T, y = T)
sel_est <- coef(sel_mod[1])[c("(Intercept)", pred_name)]
sel_est[is.na(sel_est)] <- 0
names(sel_est) <- c("(Intercept)", pred_name)
sel_se <- coef(summary(sel_mod))[, "Std. Error"][c("(Intercept)", pred_name)]
sel_se[is.na(sel_se)] <- 0

# Bootstrap results ------------------------------------------------
## The following codes were modified from Heinze G. et al (https://doi.org/10.1002/bimj.201700067)
bootnum <- 1000
boot_var <- matrix(0, ncol = length(pred) + 1, nrow = bootnum, dimnames = list(NULL, c(pred, "aic")))
boot_est <-  boot_se <- matrix(0, ncol = length(pred_name) + 1, nrow = bootnum,
                               dimnames = list(NULL, c("(Intercept)", pred_name)))
load("5a.RData") # Load bootstrap results
boot_01 <- (boot_est != 0) * 1
boot_inclusion <- apply(boot_01, 2, function(x) sum(x) / length(x) * 100)

## Overview of estimates and measures
sqe <- (t(boot_est) - full_est) ^ 2
rmsd <- apply(sqe, 1, function(x) sqrt(mean(x)))
rmsdratio <- rmsd / full_se
boot_mean <- apply(boot_est, 2, mean)
boot_meanratio <- boot_mean / full_est
boot_relbias <- (boot_meanratio / (boot_inclusion / 100) - 1) * 100
boot_median <- apply(boot_est, 2, median)
boot_025per <- apply(boot_est, 2, function(x) quantile(x, 0.025))
boot_975per <- apply(boot_est, 2, function(x) quantile(x, 0.975))
  
overview <- round(cbind(full_est, full_se, boot_inclusion, sel_est, sel_se,  
                        rmsdratio, boot_relbias, boot_median, boot_025per, 
                        boot_975per), 4)
overview <- overview[order(overview[,"boot_inclusion"], decreasing=T),]
knitr::kable(overview)
```

|             | full\_est | full\_se | boot\_inclusion | sel\_est | sel\_se | rmsdratio | boot\_relbias | boot\_median | boot\_025per | boot\_975per |
|:------------|----------:|---------:|----------------:|---------:|--------:|----------:|--------------:|-------------:|-------------:|-------------:|
| (Intercept) |  -20.6801 |   3.0499 |           100.0 | -15.7908 |  1.8360 |    1.1330 |       -1.4929 |     -20.3233 |     -26.6990 |     -14.0696 |
| IL1RA       |    0.7604 |   0.1443 |           100.0 |   0.5344 |  0.0892 |    1.2352 |        1.6338 |       0.7651 |       0.4408 |       1.1236 |
| Ang         |    0.5203 |   0.1542 |            98.2 |   0.1500 |  0.0830 |    1.1710 |        7.6528 |       0.5416 |       0.2268 |       0.8804 |
| IL8         |   -0.4250 |   0.1338 |            97.2 |  -0.0238 |  0.0851 |    1.0597 |        4.5484 |      -0.4341 |      -0.6995 |       0.0000 |
| Fer         |    0.2951 |   0.1044 |            93.5 |   0.0226 |  0.0612 |    1.2233 |       11.1579 |       0.3118 |       0.0000 |       0.5358 |
| IP10        |   -0.2145 |   0.1079 |            75.9 |  -0.1243 |  0.0586 |    1.3272 |       32.1767 |      -0.2398 |      -0.4676 |       0.0000 |
| SDC         |    0.5079 |   0.2487 |            73.5 |   0.7102 |  0.1543 |    1.2692 |       18.5555 |       0.4923 |       0.0000 |       1.0035 |
| CD163       |    0.1778 |   0.1433 |            45.8 |   0.0000 |  0.0000 |    1.1935 |       69.8797 |       0.0000 |       0.0000 |       0.4926 |
| VCAM        |   -0.0500 |   0.0541 |            35.2 |   0.0000 |  0.0000 |    1.1525 |      125.6348 |       0.0000 |      -0.1763 |       0.0000 |
| TREM        |   -0.0689 |   0.1375 |            20.0 |   0.0000 |  0.0000 |    0.9176 |      159.0696 |       0.0000 |      -0.3584 |       0.2162 |
| CRP         |    0.0438 |   0.0712 |            19.2 |   0.0000 |  0.0000 |    0.8804 |      173.2343 |       0.0000 |       0.0000 |       0.1754 |

########################################################################################### 

#### Appendix 7—table 3. Model selection frequencies for adults.

``` r
# EPV --------------------------------------------------------------
pred <- c("VCAM", "SDC", "Ang", "IL8", "ns1(IP10)", "IL1RA", "CD163", "TREM", "Fer", "CRP")

# Bootstrap results ------------------------------------------------
## The following codes were modified from Heinze G. et al (https://doi.org/10.1002/bimj.201700067)
bootnum <- 1000
boot_var <- matrix(0, ncol = length(pred) + 1, nrow = bootnum, dimnames = list(NULL, c(pred, "aic")))
load("5c.RData") # Load bootstrap results

## Model frequency
group_cols <- c(colnames(data.frame(boot_var))[1:length(pred)])
boot_modfreq <- data.frame(boot_var) %>%
  mutate(time = 1) %>%
  group_by_(.dots = group_cols) %>%
  summarise(aic_median = round(median(aic),1),
            aic_025 = round(quantile(aic, 0.025),1),
            aic_975 = round(quantile(aic, 0.975),1),
            count = sum(time)) %>%
  ungroup() %>%
  mutate(percent = count/bootnum * 100) %>%
  arrange(desc(count)) %>%
  slice(1:20) %>%
  as.data.frame(.)

out <- boot_modfreq
for (i in 1:(ncol(out)-5)) {
  out[,i] <- ifelse(out[,i]==0, NA, "+")
}

knitr::kable(out)
```

| VCAM | SDC | Ang | IL8 | ns1.IP10. | IL1RA | CD163 | TREM | Fer | CRP | aic\_median | aic\_025 | aic\_975 | count | percent |
|:-----|:----|:----|:----|:----------|:------|:------|:-----|:----|:----|------------:|---------:|---------:|------:|--------:|
|      | \+  |     | \+  | \+        | \+    | \+    | \+   | \+  |     |       407.3 |    381.3 |    442.8 |    79 |     7.9 |
| \+   | \+  |     | \+  | \+        | \+    | \+    | \+   | \+  |     |       415.0 |    374.8 |    444.2 |    55 |     5.5 |
| \+   | \+  |     | \+  | \+        | \+    |       | \+   | \+  |     |       418.0 |    380.9 |    457.7 |    36 |     3.6 |
|      | \+  |     | \+  |           | \+    |       | \+   | \+  |     |       420.5 |    389.1 |    451.4 |    33 |     3.3 |
|      | \+  |     | \+  |           |       |       | \+   | \+  | \+  |       415.6 |    371.4 |    445.2 |    30 |     3.0 |
|      | \+  | \+  | \+  | \+        | \+    | \+    | \+   | \+  |     |       406.1 |    374.3 |    434.3 |    29 |     2.9 |
|      | \+  |     | \+  |           |       | \+    | \+   | \+  | \+  |       410.3 |    375.1 |    434.6 |    26 |     2.6 |
|      | \+  |     | \+  |           | \+    |       |      | \+  |     |       428.9 |    390.0 |    460.7 |    25 |     2.5 |
|      | \+  |     | \+  |           | \+    |       | \+   | \+  | \+  |       403.8 |    388.6 |    438.9 |    24 |     2.4 |
|      | \+  |     | \+  | \+        |       | \+    | \+   | \+  | \+  |       421.0 |    384.4 |    443.8 |    20 |     2.0 |
| \+   | \+  |     | \+  |           | \+    | \+    | \+   | \+  |     |       405.8 |    384.2 |    432.4 |    20 |     2.0 |
| \+   | \+  |     | \+  | \+        | \+    | \+    | \+   | \+  | \+  |       404.9 |    374.6 |    445.0 |    20 |     2.0 |
| \+   | \+  |     | \+  | \+        | \+    | \+    |      | \+  |     |       403.1 |    384.9 |    426.5 |    19 |     1.9 |
|      | \+  |     | \+  | \+        |       | \+    | \+   | \+  |     |       424.0 |    401.4 |    450.9 |    17 |     1.7 |
|      | \+  |     | \+  |           |       |       |      | \+  | \+  |       436.1 |    409.6 |    469.2 |    16 |     1.6 |
|      | \+  |     | \+  |           |       | \+    |      | \+  |     |       422.7 |    388.3 |    455.4 |    16 |     1.6 |
|      | \+  |     | \+  |           |       | \+    | \+   | \+  |     |       416.3 |    396.7 |    465.2 |    16 |     1.6 |
|      | \+  |     | \+  | \+        |       | \+    |      | \+  |     |       417.4 |    388.5 |    451.2 |    16 |     1.6 |
|      | \+  |     | \+  | \+        | \+    |       | \+   | \+  |     |       417.9 |    396.0 |    457.2 |    16 |     1.6 |
|      | \+  |     | \+  |           | \+    | \+    | \+   | \+  |     |       412.1 |    389.1 |    435.4 |    15 |     1.5 |

########################################################################################### 

#### Appendix 7—table 4. Model stability for adults.

``` r
# EPV --------------------------------------------------------------
pred <- c("VCAM", "SDC", "Ang", "IL8", "ns1(IP10)", "IL1RA", "CD163", "TREM", "Fer", "CRP")

# Estimate full model ----------------------------------------------
full_mod <- glm(sev.or.inte ~ VCAM + SDC + Ang + IL8 + ns1(IP10) + IL1RA + CD163 + TREM + Fer + CRP, 
                family=binomial, data=dat1a, x=T, y=T)
full_est <- coef(full_mod)
full_se <- coef(summary(full_mod))[, "Std. Error"]
pred_name <- names(full_est)[-1]

# Selected model (for bootstrap results) ---------------------------
sel_m <- dredge(full_mod, rank="AIC")
sel_var <- matrix(0, ncol = length(pred), nrow = 1, dimnames = list(NULL, pred))
sel_var_tmp <- attr(get.models(sel_m, 1)[[1]]$terms, "term.labels")

for (j in 1:ncol(sel_var)) {
  sel_var[1,j] <- ifelse(names(sel_var[1,j]) %in% sel_var_tmp, 1, 0)
}

formula <- paste("sev.or.inte~", paste(names(sel_var[1,][sel_var[1,]==1]), collapse = "+"))
sel_mod <- glm(formula, data = dat1, family = binomial, x = T, y = T)
sel_est <- coef(sel_mod[1])[c("(Intercept)", pred_name)]
sel_est[is.na(sel_est)] <- 0
names(sel_est) <- c("(Intercept)", pred_name)
sel_se <- coef(summary(sel_mod))[, "Std. Error"][c("(Intercept)", pred_name)]
sel_se[is.na(sel_se)] <- 0

# Bootstrap results ------------------------------------------------
## The following codes were modified from Heinze G. et al (https://doi.org/10.1002/bimj.201700067)
bootnum <- 1000
boot_var <- matrix(0, ncol = length(pred) + 1, nrow = bootnum, dimnames = list(NULL, c(pred, "aic")))
boot_est <-  boot_se <- matrix(0, ncol = length(pred_name) + 1, nrow = bootnum,
                               dimnames = list(NULL, c("(Intercept)", pred_name)))
load("5c.RData") # Load bootstrap results
boot_01 <- (boot_est != 0) * 1
boot_inclusion <- apply(boot_01, 2, function(x) sum(x) / length(x) * 100)

## Overview of estimates and measures
sqe <- (t(boot_est) - full_est) ^ 2
rmsd <- apply(sqe, 1, function(x) sqrt(mean(x)))
rmsdratio <- rmsd / full_se
boot_mean <- apply(boot_est, 2, mean)
boot_meanratio <- boot_mean / full_est
boot_relbias <- (boot_meanratio / (boot_inclusion / 100) - 1) * 100
boot_median <- apply(boot_est, 2, median)
boot_025per <- apply(boot_est, 2, function(x) quantile(x, 0.025))
boot_975per <- apply(boot_est, 2, function(x) quantile(x, 0.975))
  
overview <- round(cbind(full_est, full_se, boot_inclusion, sel_est, sel_se,  
                        rmsdratio, boot_relbias, boot_median, boot_025per, 
                        boot_975per), 4)
overview <- overview[order(overview[,"boot_inclusion"], decreasing=T),]
knitr::kable(overview)
```

|             | full\_est | full\_se | boot\_inclusion | sel\_est | sel\_se | rmsdratio | boot\_relbias | boot\_median | boot\_025per | boot\_975per |
|:------------|----------:|---------:|----------------:|---------:|--------:|----------:|--------------:|-------------:|-------------:|-------------:|
| (Intercept) |  -16.6085 |   3.2632 |           100.0 | -16.1181 |  2.0300 |    1.3263 |        2.0586 |     -16.7262 |     -26.4755 |      -9.8095 |
| SDC         |    1.1524 |   0.2708 |            99.2 |   0.8260 |  0.1575 |    1.2745 |        3.0018 |       1.1767 |       0.5450 |       1.8615 |
| IL8         |    0.5473 |   0.1427 |            98.9 |   0.0378 |  0.0842 |    1.2973 |        5.8880 |       0.5739 |       0.2490 |       0.9520 |
| Fer         |   -0.2785 |   0.0916 |            94.6 |  -0.0083 |  0.0619 |    1.2479 |        7.9330 |      -0.2845 |      -0.5012 |       0.0000 |
| TREM        |   -0.2961 |   0.1499 |            66.5 |  -0.2171 |  0.0914 |    1.4191 |       28.1687 |      -0.2807 |      -0.6448 |       0.0000 |
| IL1RA       |    0.2582 |   0.1583 |            62.3 |   0.5217 |  0.0961 |    1.4788 |       54.9829 |       0.2494 |       0.0000 |       0.7255 |
| ns1(IP10)1  |   -1.4427 |   1.0592 |            59.8 |  -1.0181 |  0.4709 |    1.6064 |       36.8995 |      -0.2108 |      -5.2628 |       0.7003 |
| ns1(IP10)2  |   -0.1027 |   0.5763 |            59.8 |   0.1119 |  0.3226 |    1.2060 |       43.1473 |       0.0000 |      -1.7021 |       1.2589 |
| CD163       |    0.2056 |   0.1287 |            59.1 |   0.1343 |  0.0869 |    1.3333 |       51.3932 |       0.2148 |       0.0000 |       0.4989 |
| CRP         |    0.0863 |   0.0990 |            36.3 |   0.0000 |  0.0000 |    1.1203 |      129.8186 |       0.0000 |       0.0000 |       0.3048 |
| VCAM        |    0.0660 |   0.0685 |            35.4 |   0.0000 |  0.0000 |    1.3808 |       98.3111 |       0.0000 |      -0.0923 |       0.2714 |
| Ang         |   -0.0246 |   0.1193 |            21.8 |   0.0000 |  0.0000 |    1.0047 |       62.8642 |       0.0000 |      -0.2792 |       0.3008 |

########################################################################################### 

#### Appendix 8—figure 1. Results from models for severe/moderate dengue including viremia as a potential biomarker.

``` r
# Set datadist for 'lrm' function [rms]
dd <- datadist(dat1); options(datadist="dd")

# For children
dd$limits["Adjust to","Age"] <- 10

dat_m1 <- get_pred1(out="sev.or.inte", bio="VCAM", age=10, dat=dat1) %>% rename(value=VCAM) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_m2 <- get_pred1(out="sev.or.inte", bio="SDC", age=10, dat=dat1) %>% rename(value=SDC) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_m3 <- get_pred1(out="sev.or.inte", bio="Ang", age=10, dat=dat1) %>% rename(value=Ang) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_m4 <- get_pred1(out="sev.or.inte", bio="IL8", age=10, dat=dat1) %>% rename(value=IL8) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_m5 <- get_pred1(out="sev.or.inte", bio="IP10", age=10, dat=dat1) %>% rename(value=IP10) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_m6 <- get_pred1(out="sev.or.inte", bio="IL1RA", age=10, dat=dat1) %>% rename(value=IL1RA) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_m7 <- get_pred1(out="sev.or.inte", bio="CD163", age=10, dat=dat1) %>% rename(value=CD163) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_m8 <- get_pred1(out="sev.or.inte", bio="TREM", age=10, dat=dat1) %>% rename(value=TREM) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_m9 <- get_pred1(out="sev.or.inte", bio="Fer", age=10, dat=dat1) %>% rename(value=Fer) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_m10 <- get_pred1(out="sev.or.inte", bio="CRP", age=10, dat=dat1) %>% rename(value=CRP) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_m11 <- get_pred1(out="sev.or.inte", bio="Vir", age=10, dat=dat1) %>% rename(value=Vir) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))

dat_mg1 <- get_pred2v(out="sev.or.inte", bio="VCAM", age=10, dat=dat1) %>% rename(value=VCAM) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_mg2 <- get_pred2v(out="sev.or.inte", bio="SDC", age=10, dat=dat1) %>% rename(value=SDC) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_mg3 <- get_pred2v(out="sev.or.inte", bio="Ang", age=10, dat=dat1) %>% rename(value=Ang) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_mg4 <- get_pred2v(out="sev.or.inte", bio="IL8", age=10, dat=dat1) %>% rename(value=IL8) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_mg5 <- get_pred2v(out="sev.or.inte", bio="IP10", age=10, dat=dat1) %>% rename(value=IP10) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_mg6 <- get_pred2v(out="sev.or.inte", bio="IL1RA", age=10, dat=dat1) %>% rename(value=IL1RA) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_mg7 <- get_pred2v(out="sev.or.inte", bio="CD163", age=10, dat=dat1) %>% rename(value=CD163) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_mg8 <- get_pred2v(out="sev.or.inte", bio="TREM", age=10, dat=dat1) %>% rename(value=TREM) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_mg9 <- get_pred2v(out="sev.or.inte", bio="Fer", age=10, dat=dat1) %>% rename(value=Fer) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_mg10 <- get_pred2v(out="sev.or.inte", bio="CRP", age=10, dat=dat1) %>% rename(value=CRP) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))
dat_mg11 <- get_pred2v(out="sev.or.inte", bio="Vir", age=10, dat=dat1) %>% rename(value=Vir) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="Under 15"))

dat_p1 <- rbind(dat_m1,dat_mg1, dat_m2,dat_mg2, dat_m3,dat_mg3, dat_m4,dat_mg4, dat_m5,dat_mg5, dat_m6,dat_mg6, dat_m7,dat_mg7, dat_m8,dat_mg8, dat_m9,dat_mg9, dat_m10,dat_mg10, dat_m11,dat_mg11) %>%
  mutate(age = "10 years")

# For adults
dd$limits["Adjust to","Age"] <- 25

dat_m1 <- get_pred1(out="sev.or.inte", bio="VCAM", age=25, dat=dat1) %>% rename(value=VCAM) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_m2 <- get_pred1(out="sev.or.inte", bio="SDC", age=25, dat=dat1) %>% rename(value=SDC) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_m3 <- get_pred1(out="sev.or.inte", bio="Ang", age=25, dat=dat1) %>% rename(value=Ang) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_m4 <- get_pred1(out="sev.or.inte", bio="IL8", age=25, dat=dat1) %>% rename(value=IL8) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_m5 <- get_pred1(out="sev.or.inte", bio="IP10", age=25, dat=dat1) %>% rename(value=IP10) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_m6 <- get_pred1(out="sev.or.inte", bio="IL1RA", age=25, dat=dat1) %>% rename(value=IL1RA) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_m7 <- get_pred1(out="sev.or.inte", bio="CD163", age=25, dat=dat1) %>% rename(value=CD163) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_m8 <- get_pred1(out="sev.or.inte", bio="TREM", age=25, dat=dat1) %>% rename(value=TREM) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_m9 <- get_pred1(out="sev.or.inte", bio="Fer", age=25, dat=dat1) %>% rename(value=Fer) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_m10 <- get_pred1(out="sev.or.inte", bio="CRP", age=25, dat=dat1) %>% rename(value=CRP) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_m11 <- get_pred1(out="sev.or.inte", bio="Vir", age=25, dat=dat1) %>% rename(value=Vir) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))

dat_mg1 <- get_pred2v(out="sev.or.inte", bio="VCAM", age=25, dat=dat1) %>% rename(value=VCAM) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_mg2 <- get_pred2v(out="sev.or.inte", bio="SDC", age=25, dat=dat1) %>% rename(value=SDC) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_mg3 <- get_pred2v(out="sev.or.inte", bio="Ang", age=25, dat=dat1) %>% rename(value=Ang) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_mg4 <- get_pred2v(out="sev.or.inte", bio="IL8", age=25, dat=dat1) %>% rename(value=IL8) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_mg5 <- get_pred2v(out="sev.or.inte", bio="IP10", age=25, dat=dat1) %>% rename(value=IP10) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_mg6 <- get_pred2v(out="sev.or.inte", bio="IL1RA", age=25, dat=dat1) %>% rename(value=IL1RA) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_mg7 <- get_pred2v(out="sev.or.inte", bio="CD163", age=25, dat=dat1) %>% rename(value=CD163) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_mg8 <- get_pred2v(out="sev.or.inte", bio="TREM", age=25, dat=dat1) %>% rename(value=TREM) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_mg9 <- get_pred2v(out="sev.or.inte", bio="Fer", age=25, dat=dat1) %>% rename(value=Fer) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_mg10 <- get_pred2v(out="sev.or.inte", bio="CRP", age=25, dat=dat1) %>% rename(value=CRP) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))
dat_mg11 <- get_pred2v(out="sev.or.inte", bio="Vir", age=25, dat=dat1) %>% rename(value=Vir) %>% 
  select(value, yhat, lower, upper, biomarker, model) %>% slice(which(dat1$age15=="15 and above"))

dat_p2 <- rbind(dat_m1,dat_mg1, dat_m2,dat_mg2, dat_m3,dat_mg3, dat_m4,dat_mg4, dat_m5,dat_mg5, dat_m6,dat_mg6, dat_m7,dat_mg7, dat_m8,dat_mg8, dat_m9,dat_mg9, dat_m10,dat_mg10, dat_m11,dat_mg11) %>%
  mutate(age = "25 years")

# Merge data for children and adults
vis1 <- rbind(dat_p1, dat_p2) %>% mutate(outcome = "outcome 1")

# Merge data for plots
tmp0 <- dat %>% arrange(age15) %>% select(sev.or.inte)
tmp <- data.frame(sev.or.inte = rep(tmp0$sev.or.inte, 22))
  
tmp_p1 <- vis1 %>%
  arrange(biomarker, model) %>%
  mutate(value1 = 2^value,
         age = factor(age, levels=c("10 years", "25 years")),
         gr = ifelse(model=="Single model" & age=="10 years", 1,
                     ifelse(model=="Single model" & age=="25 years", 2,
                            ifelse(model=="Global model" & age=="10 years", 3, 4))),
         gr = factor(gr, levels=c(1:4), labels=c("Single children", "Single adults", "Global children", "Global adults")),
         model = factor(model, levels=c("Single model", "Global model"))) %>%
  bind_cols(., tmp)

vline <- tmp_p1 %>%
  filter(model=="Single model") %>%
  filter(yhat==0) %>%
  select(biomarker, value, value1) %>%
  rename(vline=value, vline1=value1)

dat_p <- left_join(tmp_p1, vline, by="biomarker")

# Alternative data to limit to 5th - 95th of the values
dat_alt <- dat_p %>%
  filter(value != 0) %>%
  group_by(biomarker, model) %>%
  mutate(up = quantile(value, .95),
         lo = quantile(value, .05)) %>%
  ungroup() %>%
  mutate(is.outlier = value<lo | value>up | value<0,
         value = ifelse(is.outlier, NA, value),
         value1 = ifelse(is.outlier, NA, value1),
         yhat = ifelse(is.outlier, NA, yhat),
         lower = ifelse(is.outlier, NA, lower),
         upper = ifelse(is.outlier, NA, upper))

dat_1_5 <- dat_alt %>% 
  filter(biomarker %in% c("VCAM", "SDC", "Ang", "IL8", "IP10")) %>%
  mutate(biomarker = factor(biomarker, levels=c("VCAM", "SDC", "Ang", "IL8", "IP10")))

dat_6_11 <- dat_alt %>% 
  filter(biomarker %in% c("IL1RA", "CD163", "TREM", "Fer", "CRP", "Vir")) %>%
  mutate(biomarker = factor(biomarker, levels=c("IL1RA", "CD163", "TREM", "Fer", "CRP", "Vir")))

# Modify facets' scales
#library(facetscales)

xVCAM <- c(1,4,15,60,250,1000,4000)
xSDC <- c(1400,2000,2800,4000,5600)
xAng <- c(50,100,250,500,1000,2000)
xIL8 <- c(5,7,10,14,20,28,40)
xIP10 <- c(25,100,400,1600,6400)
xIL1RA <- c(1000,2000,4000,8000,16000)
xCD163 <- c(75,150,300,600)
xTREM <- c(35,50,70,100,140,200)
xFer <- c(50,100,200,400,800)
xCRP <- c(2.5,5,10,20,40,80)
xVir <- c(0:10)

scales_x <- list(
  `VCAM` = scale_x_continuous(breaks = log2(xVCAM), labels = xVCAM),
  `SDC` = scale_x_continuous(breaks = log2(xSDC), labels = xSDC),
  `Ang` = scale_x_continuous(breaks = log2(xAng), labels = xAng), 
  `IL8` = scale_x_continuous(breaks = log2(xIL8), labels = xIL8), 
  `IP10` = scale_x_continuous(breaks = log2(xIP10), labels = xIP10), 
  `IL1RA` = scale_x_continuous(breaks = log2(xIL1RA), labels = xIL1RA), 
  `CD163` = scale_x_continuous(breaks = log2(xCD163), labels = xCD163), 
  `TREM` = scale_x_continuous(breaks = log2(xTREM), labels = xTREM), 
  `Fer` = scale_x_continuous(breaks = log2(xFer), labels = xFer), 
  `CRP` = scale_x_continuous(breaks = log2(xCRP), labels = xCRP),
  `Vir` = scale_x_continuous(breaks = xVir, labels = xVir)
)

ybreak1 <- c(.125, .25, .5, 1, 2, 4, 8)
ybreak2 <- c(.125, .25, .5, 1, 2, 4)

scales_y1 <- list(
  `Single model` = scale_y_continuous(limits = c(log(.124), log(8.01)), breaks = log(ybreak1), labels = ybreak1),
  `Global model` = scale_y_continuous(limits = c(log(.124), log(8.01)), breaks = log(ybreak1), labels = ybreak1)
)

scales_y2 <- list(
  `Single model` = scale_y_continuous(limits = c(log(.124), log(4.01)), breaks = log(ybreak2), labels = ybreak2),
  `Global model` = scale_y_continuous(limits = c(log(.124), log(4.01)), breaks = log(ybreak2), labels = ybreak2)
)

# Set facets' names
models <- c(`Single model` = "Single model",
            `Global model` = "Global model")

biomarkers <- c(`VCAM` = "VCAM-1 (ng/ml)",
                `SDC`= "SDC-1 (pg/ml)",
                `Ang` = "Ang-2 (pg/ml)",
                `IL8` = "IL-8 (pg/ml)",
                `IP10` = "IP-10 (pg/ml)",
                `IL1RA` = "IL-1RA (pg/ml)",
                `CD163` = "sCD163 (ng/ml)",
                `TREM` = "sTREM-1 (pg/ml)",
                `Fer` = "Ferritin (ng/ml)",
                `CRP` = "CRP (mg/l)",
                `Vir` = "Log10 viremia")

# Plot for the first 5 biomarkers
p.1.5 <- ggplot(dat_1_5, aes(x=value)) +
  geom_ribbon(aes(ymin=lower, ymax=upper, fill=age), alpha=.15) +
  geom_vline(aes(xintercept=vline), color="black", linetype="dashed", alpha=.4) +
  geom_hline(yintercept=0, color="black", linetype="dashed", alpha=.4) +
  geom_line(aes(y=yhat, color=age), size=.5, alpha=.7) +
  geom_rug(data=filter(dat_1_5, model=="Single model" & sev.or.inte==0 & age=="10 years"), sides="b", alpha=.2, color="red") +
  geom_rug(data=filter(dat_1_5, model=="Single model" & sev.or.inte==1 & age=="10 years"), sides="t", alpha=.2, color="red") +
  geom_rug(data=filter(dat_1_5, model=="Global model" & sev.or.inte==0 & age=="25 years"), sides="b", alpha=.2, color="blue") +
  geom_rug(data=filter(dat_1_5, model=="Global model" & sev.or.inte==1 & age=="25 years"), sides="t", alpha=.2, color="blue") +
  scale_color_manual(values=c("red","blue"), labels=c("Children","Adults")) +
  scale_fill_manual(values=c("red","blue"), labels=c("Children","Adults")) +
  ylab("Odds ratio") +
  theme(legend.position="top", legend.title=element_blank(), axis.title.x=element_blank(), axis.text.x=element_text(angle=45)) +
  facet_grid_sc(cols=vars(biomarker), rows=vars(model), 
                scales=list(x=scales_x, y=scales_y1),
                labeller = as_labeller(c(models, biomarkers)))

# Plot for the last 5 biomarkers
p.6.11 <- ggplot(dat_6_11, aes(x=value)) +
  geom_ribbon(aes(ymin=lower, ymax=upper, fill=age), alpha=.15) +
  geom_vline(aes(xintercept=vline), color="black", linetype="dashed", alpha=.4) +
  geom_hline(yintercept=0, color="black", linetype="dashed", alpha=.4) +
  geom_line(aes(y=yhat, color=age), size=.5, alpha=.7) +
  geom_rug(data=filter(dat_6_11, model=="Single model" & sev.or.inte==0 & age=="10 years"), sides="b", alpha=.2, color="red") +
  geom_rug(data=filter(dat_6_11, model=="Single model" & sev.or.inte==1 & age=="10 years"), sides="t", alpha=.2, color="red") +
  geom_rug(data=filter(dat_6_11, model=="Global model" & sev.or.inte==0 & age=="25 years"), sides="b", alpha=.2, color="blue") +
  geom_rug(data=filter(dat_6_11, model=="Global model" & sev.or.inte==1 & age=="25 years"), sides="t", alpha=.2, color="blue") +
  scale_color_manual(values=c("red","blue"), labels=c("Children","Adults")) +
  scale_fill_manual(values=c("red","blue"), labels=c("Children","Adults")) +
  ylab("Odds ratio") +
  theme(legend.position="none", axis.title.x=element_blank(), axis.text.x=element_text(angle=45)) +
  facet_grid_sc(cols=vars(biomarker), rows=vars(model), 
                scales=list(x=scales_x, y=scales_y2),
                labeller = as_labeller(c(models, biomarkers)))

# Merge plots
p81 <- gridExtra::grid.arrange(p.1.5, p.6.11, nrow=2, heights=c(1,.86))
```

![](Elife-ERA-codes_files/figure-gfm/fig%20A8.1-1.png)<!-- -->

``` r
#ggsave(filename="A8_fig1.pdf", p81, dpi=300, width=8.5, height=6.2)
```

########################################################################################### 

#### Appendix 8—table 1. Results from models for severe/moderate dengue including viremia as a potential biomarker.

``` r
# Set datadist for 'lrm' function [rms]
dd <- datadist(dat1); options(datadist="dd")

# Set references to calculate results
ref0v <- c(ref0, median(dat1$Vir))

# Get ORs and 95% CIs
tmp1 <- data.frame(
  sort1 = rep(c(1:11), 2),
  sort2 = c(rep(2,11), rep(3,11)),
  bio = rep(c("VCAM", "SDC", "Ang", "IL8", "IP10", "IL1RA", "CD163", "TREM", "Fer", "CRP", "Vir"), 2),
  ref1 = c(ref0v-1, ref0v),
  ref2 = c(ref0v, ref0v+1)
) %>%
  arrange(sort1, sort2) %>%
  group_by(sort1, sort2) %>%
  do(cbind(.,
           # Children - single models
           or1c = get_est1(out="sev.or.inte", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="OR", age=10, dat=dat1),
           lo1c = get_est1(out="sev.or.inte", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="loCI", age=10, dat=dat1),
           up1c = get_est1(out="sev.or.inte", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="upCI", age=10, dat=dat1),
           # Children - global model
           or2c = get_est2v(out="sev.or.inte", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="OR", age=10, dat=dat1),
           lo2c = get_est2v(out="sev.or.inte", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="loCI", age=10, dat=dat1),
           up2c = get_est2v(out="sev.or.inte", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="upCI", age=10, dat=dat1),
           # Adults - single models
           or1a = get_est1(out="sev.or.inte", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="OR", age=25, dat=dat1),
           lo1a = get_est1(out="sev.or.inte", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="loCI", age=25, dat=dat1),
           up1a = get_est1(out="sev.or.inte", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="upCI", age=25, dat=dat1),
           # Adults - global model
           or2a = get_est2v(out="sev.or.inte", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="OR", age=25, dat=dat1),
           lo2a = get_est2v(out="sev.or.inte", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="loCI", age=25, dat=dat1),
           up2a = get_est2v(out="sev.or.inte", bio=.$bio, ref1=.$ref1, ref2=.$ref2, est="upCI", age=25, dat=dat1))) %>%
  ungroup()
for (i in 6:17) {tmp1[[i]] <- sprintf("%.2f", round(tmp1[[i]],2))}

# Get p-values
tmp2 <- data.frame(
  sort1 = c(1:11),
  sort2 = rep(1,11),
  bio = c("VCAM", "SDC", "Ang", "IL8", "IP10", "IL1RA", "CD163", "TREM", "Fer", "CRP", "Vir")
) %>%
  group_by(sort1) %>%
  do(cbind(.,
           p.s1 = get_est1(out="sev.or.inte", bio=.$bio, est="p", dat=dat1),
           p.s1_int = get_est1(out="sev.or.inte", bio=.$bio, est="p int", dat=dat1),
           p.g1 = get_est2v(out="sev.or.inte", bio=.$bio, est="p", dat=dat1),
           p.g1_int = get_est2v(out="sev.or.inte", bio=.$bio, est="p int", dat=dat1))) %>%
  ungroup()
for (i in 4:7) {tmp2[[i]] <- ifelse(tmp2[[i]]<0.001, "<0.001", sprintf("%.3f", round(tmp2[[i]],3)))}  

# Combine results
res1 <- bind_rows(tmp1, tmp2) %>%
  arrange(sort1, sort2) %>%
  mutate(bio = ifelse(is.na(ref1), as.character(bio),
                      ifelse(bio!="Vir", paste(" -", round(2^ref2,0), "vs", round(2^ref1,0), sep=" "),
                             paste(" -", round(ref2,1), "vs", round(ref1,1), sep=" "))),
         or.sc1 = ifelse(is.na(lo1c), NA, paste(or1c, " (", lo1c, "-", up1c, ")", sep="")), # s: single model; c: children
         or.gc1 = ifelse(is.na(lo2c), NA, paste(or2c, " (", lo2c, "-", up2c, ")", sep="")), # g: global model
         or.sa1 = ifelse(is.na(lo1a), NA, paste(or1a, " (", lo1a, "-", up1a, ")", sep="")), # a: adults
         or.ga1 = ifelse(is.na(lo2a), NA, paste(or2a, " (", lo2a, "-", up2a, ")", sep=""))) %>%
  select(bio, or.sc1, or.sa1, p.s1, p.s1_int, or.gc1, or.ga1, p.g1, p.g1_int)

names(res1) <- c("", "OR (single - children)", "OR (single - adults)", "P overall (single)", "P interaction (single)",
                 "OR (global - children)", "OR (global - adults)", "P overall (global)", "P interaction (global)")

# Report results
knitr::kable(res1)
```

|                  | OR (single - children) | OR (single - adults) | P overall (single) | P interaction (single) | OR (global - children) | OR (global - adults) | P overall (global) | P interaction (global) |
|:-----------------|:-----------------------|:---------------------|:-------------------|:-----------------------|:-----------------------|:---------------------|:-------------------|:-----------------------|
| VCAM             |                        |                      | &lt;0.001          | 0.715                  |                        |                      | 0.286              | 0.136                  |
| \- 1636 vs 818   | 1.20 (1.04-1.38)       | 1.35 (1.15-1.58)     |                    |                        | 0.94 (0.76-1.17)       | 1.34 (1.03-1.74)     |                    |                        |
| \- 3272 vs 1636  | 1.25 (1.02-1.53)       | 1.48 (1.19-1.85)     |                    |                        | 0.93 (0.69-1.24)       | 1.45 (1.02-2.04)     |                    |                        |
| SDC              |                        |                      | &lt;0.001          | 0.088                  |                        |                      | 0.005              | 0.645                  |
| \- 2519 vs 1260  | 2.67 (1.31-5.43)       | 3.33 (1.32-8.42)     |                    |                        | 2.07 (0.78-5.47)       | 4.28 (1.27-14.43)    |                    |                        |
| \- 5039 vs 2519  | 1.71 (1.18-2.47)       | 3.71 (2.09-6.58)     |                    |                        | 1.71 (0.95-3.09)       | 2.55 (1.17-5.57)     |                    |                        |
| Ang              |                        |                      | &lt;0.001          | 0.524                  |                        |                      | 0.060              | 0.070                  |
| \- 1204 vs 602   | 1.64 (1.39-1.94)       | 1.51 (1.26-1.82)     |                    |                        | 1.62 (1.19-2.20)       | 0.97 (0.71-1.34)     |                    |                        |
| \- 2409 vs 1204  | 2.21 (1.58-3.10)       | 2.00 (1.40-2.85)     |                    |                        | 1.92 (1.22-3.01)       | 0.96 (0.61-1.49)     |                    |                        |
| IL8              |                        |                      | &lt;0.001          | &lt;0.001              |                        |                      | &lt;0.001          | &lt;0.001              |
| \- 14 vs 7       | 1.42 (1.05-1.91)       | 2.18 (1.47-3.24)     |                    |                        | 0.89 (0.61-1.31)       | 1.60 (0.99-2.59)     |                    |                        |
| \- 28 vs 14      | 0.99 (0.78-1.25)       | 2.33 (1.63-3.33)     |                    |                        | 0.52 (0.36-0.77)       | 1.96 (1.28-3.02)     |                    |                        |
| IP10             |                        |                      | &lt;0.001          | 0.984                  |                        |                      | 0.150              | 0.500                  |
| \- 3093 vs 1546  | 1.46 (1.26-1.68)       | 1.45 (1.21-1.73)     |                    |                        | 0.93 (0.73-1.19)       | 0.75 (0.53-1.06)     |                    |                        |
| \- 6186 vs 3093  | 1.68 (1.35-2.09)       | 1.69 (1.29-2.22)     |                    |                        | 1.07 (0.76-1.50)       | 0.75 (0.47-1.20)     |                    |                        |
| IL1RA            |                        |                      | &lt;0.001          | 0.082                  |                        |                      | &lt;0.001          | 0.062                  |
| \- 6434 vs 3217  | 1.69 (1.42-2.03)       | 1.48 (1.21-1.81)     |                    |                        | 1.97 (1.42-2.73)       | 1.40 (0.93-2.09)     |                    |                        |
| \- 12868 vs 6434 | 1.82 (1.46-2.27)       | 1.70 (1.29-2.24)     |                    |                        | 2.03 (1.41-2.92)       | 1.38 (0.87-2.19)     |                    |                        |
| CD163            |                        |                      | &lt;0.001          | 0.551                  |                        |                      | 0.124              | 0.289                  |
| \- 295 vs 147    | 1.57 (1.14-2.15)       | 1.49 (1.13-1.98)     |                    |                        | 1.51 (0.94-2.42)       | 1.30 (0.86-1.98)     |                    |                        |
| \- 589 vs 295    | 1.46 (1.10-1.93)       | 1.61 (1.09-2.37)     |                    |                        | 1.24 (0.88-1.73)       | 1.44 (0.91-2.28)     |                    |                        |
| TREM             |                        |                      | 0.059              | 0.997                  |                        |                      | 0.745              | 0.594                  |
| \- 85 vs 42      | 1.87 (1.23-2.84)       | 1.79 (1.10-2.93)     |                    |                        | 1.16 (0.71-1.91)       | 1.24 (0.65-2.36)     |                    |                        |
| \- 169 vs 85     | 1.12 (0.91-1.38)       | 1.12 (0.82-1.53)     |                    |                        | 0.92 (0.67-1.27)       | 0.67 (0.41-1.10)     |                    |                        |
| Fer              |                        |                      | 0.042              | 0.054                  |                        |                      | 0.007              | 0.002                  |
| \- 243 vs 122    | 1.18 (1.01-1.38)       | 1.06 (0.89-1.27)     |                    |                        | 1.32 (1.04-1.67)       | 0.77 (0.60-0.98)     |                    |                        |
| \- 487 vs 243    | 1.26 (1.00-1.58)       | 0.90 (0.66-1.23)     |                    |                        | 1.22 (0.89-1.68)       | 0.64 (0.42-0.98)     |                    |                        |
| CRP              |                        |                      | &lt;0.001          | 0.031                  |                        |                      | 0.185              | 0.113                  |
| \- 28 vs 14      | 1.26 (1.12-1.41)       | 1.25 (1.03-1.52)     |                    |                        | 1.06 (0.91-1.23)       | 1.09 (0.84-1.43)     |                    |                        |
| \- 56 vs 28      | 1.13 (0.95-1.34)       | 1.38 (1.11-1.71)     |                    |                        | 0.89 (0.72-1.11)       | 1.35 (1.01-1.81)     |                    |                        |
| Vir              |                        |                      | &lt;0.001          | 0.747                  |                        |                      | 0.040              | 0.886                  |
| \- 7.5 vs 6.5    | 1.34 (1.18-1.53)       | 1.33 (1.16-1.53)     |                    |                        | 1.21 (1.03-1.42)       | 1.25 (1.04-1.51)     |                    |                        |
| \- 8.5 vs 7.5    | 1.53 (1.25-1.87)       | 1.48 (1.16-1.88)     |                    |                        | 1.35 (1.05-1.74)       | 1.43 (1.05-1.95)     |                    |                        |

########################################################################################### 

#### Appendix 8—table 2. Best combinations of biomarkers associated with severe or moderate dengue for children.

``` r
# EPV --------------------------------------------------------------
pred <- c("VCAM", "SDC", "Ang", "IL8", "IP10", "IL1RA", "CD163", "TREM", "Fer", "CRP", "Vir")

# Estimate full model ----------------------------------------------
full_mod <- glm(sev.or.inte ~ VCAM + SDC + Ang + IL8 + IP10 + IL1RA + CD163 + TREM + Fer + CRP + Vir, family=binomial, data=dat1c, x=T, y=T)

# Selected model ---------------------------------------------------
sel_var <- matrix(0, ncol=length(pred)+1, nrow=5, dimnames=list(NULL, c(pred, "aic")))

for (i in 1:5) {
  if (i==1) {
    bs <- dredge(full_mod, rank="AIC")
  } else {
    bs <- dredge(full_mod, rank="AIC", m.lim=c(i,i))
  }
  bs_var <- attr(get.models(bs, 1)[[1]]$terms, "term.labels")

  for (j in 1:(ncol(sel_var)-1)) {sel_var[i,j] <- ifelse(names(sel_var[i,j]) %in% bs_var, 1, 0)}
  formula <- paste("sev.or.inte~", paste(names(sel_var[i,][sel_var[i,]==1]), collapse = "+"))
  sel_mod <- glm(formula, data = dat1c, family = binomial, x = T, y = T)
  sel_var[i, ncol(sel_var)] <- AIC(sel_mod)
}

# Report results --------------------------------------------------
out1 <- as.data.frame(sel_var) %>% mutate(AIC = round(aic,1)) %>% select(-aic)
for (i in 1:(ncol(out1)-1)) {out1[,i] <- ifelse(out1[,i]==0, NA, "+")}

out2 <- as.data.frame(t(out1))
colnames(out2) <- c("Best of all combinations", "Best combination of 2 variables", "Best combination of 3 variables",
                    "Best combination of 4 variables", "Best combination of 5 variables")
rownames(out2) <- c("- VCAM-1", "- SDC-1", "- Ang-2", "- IL-8", "- IP-10", "- IL-1RA", "- sCD163", "- sTREM-1",
                    "- Ferritin", "- CRP", "- Viremia", "AIC of the selected model")

knitr::kable(out2)
```

|                           | Best of all combinations | Best combination of 2 variables | Best combination of 3 variables | Best combination of 4 variables | Best combination of 5 variables |
|:--------------------------|:-------------------------|:--------------------------------|:--------------------------------|:--------------------------------|:--------------------------------|
| \- VCAM-1                 |                          |                                 |                                 |                                 |                                 |
| \- SDC-1                  | \+                       |                                 |                                 |                                 |                                 |
| \- Ang-2                  | \+                       |                                 | \+                              | \+                              | \+                              |
| \- IL-8                   | \+                       |                                 |                                 |                                 | \+                              |
| \- IP-10                  | \+                       |                                 |                                 | \+                              | \+                              |
| \- IL-1RA                 | \+                       | \+                              | \+                              | \+                              | \+                              |
| \- sCD163                 |                          |                                 |                                 |                                 |                                 |
| \- sTREM-1                |                          |                                 |                                 |                                 |                                 |
| \- Ferritin               | \+                       | \+                              | \+                              | \+                              | \+                              |
| \- CRP                    |                          |                                 |                                 |                                 |                                 |
| \- Viremia                |                          |                                 |                                 |                                 |                                 |
| AIC of the selected model | 465.9                    | 484.7                           | 480.0                           | 473.7                           | 467.6                           |

########################################################################################### 

#### Appendix 8—table 3. Best combinations of biomarkers associated with severe or moderate dengue for adults.

``` r
# EPV --------------------------------------------------------------
pred <- c("VCAM", "SDC", "Ang", "IL8", "ns1(IP10)", "IL1RA", "CD163", "TREM", "Fer", "CRP", "Vir")

# Estimate full model ----------------------------------------------
full_mod <- glm(sev.or.inte ~ VCAM + SDC + Ang + IL8 + ns1(IP10) + IL1RA + CD163 + TREM + Fer + CRP + Vir, family=binomial, data=dat1a, x=T, y=T)

# Selected model ---------------------------------------------------
sel_var <- matrix(0, ncol=length(pred)+1, nrow=5, dimnames=list(NULL, c(pred, "aic")))

for (i in 1:5) {
  if (i==1) {
    bs <- dredge(full_mod, rank="AIC")
  } else {
    bs <- dredge(full_mod, rank="AIC", m.lim=c(i,i))
  }
  bs_var <- attr(get.models(bs, 1)[[1]]$terms, "term.labels")

  for (j in 1:(ncol(sel_var)-1)) {sel_var[i,j] <- ifelse(names(sel_var[i,j]) %in% bs_var, 1, 0)}
  formula <- paste("sev.or.inte~", paste(names(sel_var[i,][sel_var[i,]==1]), collapse = "+"))
  sel_mod <- glm(formula, data = dat1a, family = binomial, x = T, y = T)
  sel_var[i, ncol(sel_var)] <- AIC(sel_mod)
}

# Report results --------------------------------------------------
out1 <- as.data.frame(sel_var) %>% mutate(AIC = round(aic,1)) %>% select(-aic)
for (i in 1:(ncol(out1)-1)) {out1[,i] <- ifelse(out1[,i]==0, NA, "+")}

out2 <- as.data.frame(t(out1))
colnames(out2) <- c("Best of all combinations", "Best combination of 2 variables", "Best combination of 3 variables",
                    "Best combination of 4 variables", "Best combination of 5 variables")
rownames(out2) <- c("- VCAM-1", "- SDC-1", "- Ang-2", "- IL-8", "- IP-10", "- IL-1RA", "- sCD163", "- sTREM-1",
                    "- Ferritin", "- CRP", "- Viremia", "AIC of the selected model")

knitr::kable(out2)
```

|                           | Best of all combinations | Best combination of 2 variables | Best combination of 3 variables | Best combination of 4 variables | Best combination of 5 variables |
|:--------------------------|:-------------------------|:--------------------------------|:--------------------------------|:--------------------------------|:--------------------------------|
| \- VCAM-1                 |                          |                                 |                                 |                                 |                                 |
| \- SDC-1                  | \+                       | \+                              | \+                              | \+                              | \+                              |
| \- Ang-2                  |                          |                                 |                                 |                                 |                                 |
| \- IL-8                   | \+                       | \+                              | \+                              | \+                              | \+                              |
| \- IP-10                  |                          |                                 |                                 |                                 |                                 |
| \- IL-1RA                 |                          |                                 |                                 |                                 |                                 |
| \- sCD163                 | \+                       |                                 |                                 |                                 | \+                              |
| \- sTREM-1                |                          |                                 |                                 |                                 |                                 |
| \- Ferritin               | \+                       |                                 | \+                              | \+                              | \+                              |
| \- CRP                    |                          |                                 |                                 |                                 |                                 |
| \- Viremia                | \+                       |                                 |                                 | \+                              | \+                              |
| AIC of the selected model | 426.4                    | 441.1                           | 434.2                           | 428.5                           | 426.4                           |

########################################################################################### 

#### Appendix 9—table 1. Results of variable selection for children.

########################################################################################### 

#### Appendix 9—table 2. Results of variable selection for adults.
