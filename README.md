# South African SARS-CoV-2 Household Transmission Study (SA-S-HTS) Proximity Analysis

This repository houses data used for Kleynhans. et al. 2023. Association of close-range contact patterns with SARS-CoV-2 household transmission, South Africa, 2020-2021.



Created by: Jackie Kleynhans, National Institute for Communicable Diseases

Households are an important location for severe acute respiratory syndrome coronavirus 2 (SARS-CoV-2) transmission, especially during periods where travel and work was restricted to essential services. We aimed to assess the association of close-range contact patterns with SARS-CoV-2 transmission.

We deployed proximity sensors for two weeks to measure face-to-face interactions between household members after SARS-CoV-2 was identified in the household, in South Africa, 2020 - 2021. We calculated duration, frequency and average duration of close range proximity events with SARS-CoV-2 index cases. We assessed the association of contact parameters with SARS-CoV-2 transmission using mixed effects logistic regression accounting for index and household member characteristics.

Details of main household transmisison study published: Kleynhans J, Walaza S, Martinson NA, Neti M, von Gottberg A, Bhiman JN, et al. 2022. Household transmission of SARS-CoV-2 from adult index cases living with and without HIV in South Africa, 2020-2021: A case-ascertained, prospective observational household transmission study. Clinical Infectious Diseases. DOI: https://doi.org/10.1093/cid/ciac640.

To start, download and unzip files from sashts_data.zip which contains two files: sashts_contact_network.csv the complete network file for household close proximity events (Symmetrical network where all close proximity events appear twice. Therefore all contacts for each individual in indid1 exists when considering only indid1 column) and sashts_metadata.csv which includes demographics of household members
  
 Data dictionary also provided: Data dictionary.xlsx is stored in main repository
  
 Analyses codes: SA-S-HTS Proximity Analysis.R run on R-4.1.1 in R Studio 2021.09.0 Build 351
