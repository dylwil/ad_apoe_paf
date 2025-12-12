

* Creating an indicator for randomly excluding related individuals from within pairs:
clear
import delimited "ukb_rel.dat", delimiter(whitespace)

 gen rand = runiform()
 
 gen exclusions = id1 if rand<=.5
 replace exclusions = id2 if rand>.5
 
 keep exclusions
 
 duplicates drop
 
 rename exclusions eid

 gen simple_rel_exclusions = 1
 
 save "rel_exclusions_simple.dta", replace
 
 clear
 
* Mark-out of all-cause dementia from primary care records (based on record extract conducted with ukbrapR package) 

use "dem_gp_records.dta"
sort eid

gen dem_pc = 1 
duplicates drop eid, force

keep eid dem_pc

save "dem_pc_markout.dta", replace


* Mark-out of AD from primary care records (based on record extract conducted with ukbrapR package) 
clear
use "alz_gp_records.dta"
sort eid

gen alz_pc = 1 
duplicates drop eid, force

keep eid alz_pc

save "alz_pc_markout.dta", replace

clear

* Importing phenotype data curated with cohort builder:
clear
import delimited "apoe_paf_pheno.csv"

gen age_baseline = p21003 if p21003>=60
keep if age_baseline!=.



* APOE (based on SNP extracts conducted with ukbrapR package) 

merge 1:1 eid using "imputed_apoe.dta" 

keep if _m==3

gen e3_e3 = 1 if rs429358_C==0 & rs7412_T==0 // TT for rs429358 and CC for rs7412
gen e2_e3 = 1 if rs429358_C==0 & rs7412_T==1 // TT for rs429358 and CT for rs7412
gen e2_e2 = 1 if rs429358_C==0 & rs7412_T==2 // TT for rs429358 and TT for rs7412

gen e3_e4 = 1 if rs429358_C==1 & rs7412_T==0 // CT for rs429358 and CC for rs7412
gen e2_e4 = 1 if rs429358_C==1 & rs7412_T==1 // CT for rs429358 and CT for rs7412 - nb: e2/e4 and e1/e3 genotypes ambiguous without phased haplotypes
gen  e1_e2 = 1 if rs429358_C==1 & rs7412_T==2 // CT for rs429358 and TT for rs7412

gen e4_e4 = 1 if rs429358_C==2 & rs7412_T==0 // CC for rs429358 and CC for rs7412
gen e1_e4 = 1 if rs429358_C==2 & rs7412_T==1 // CC for rs429358 and CT for rs7412
gen  e1_e1 = 1 if rs429358_C==2 & rs7412_T==2 // CC for rs429358 and TT for rs7412

foreach var in e1_e4 e1_e2 e2_e4 e2_e2 e3_e3 e3_e4 e4_e4 e2_e3 {
tab `var'
}

gen apoe_genotype = 1 if e1_e4==1
replace apoe_genotype = 2 if e1_e2==1
replace apoe_genotype = 3 if e2_e2==1
replace apoe_genotype = 4 if e2_e3==1
replace apoe_genotype = 5 if e3_e3==1 
replace apoe_genotype = 6 if e2_e4==1
replace apoe_genotype = 7 if e3_e4==1
replace apoe_genotype = 8 if e4_e4==1

label define apoe 1 "e1_e4" 2 "e1_e2" 3 "e2_e2" 4 "e2_e3" 5 "e3_e3" 6 "e2_e4" 7 "e3_e4" 8 "e4_e4"
label values apoe_genotype  apoe

* Removal of individuals carrying presumed e1 genotypes
drop if apoe_genotype<3

gen e3_carrier = 0
replace e3_carrier = 1 if e2_e3==1
replace e3_carrier = 1 if e3_e3==1
replace e3_carrier = 1 if e3_e4==1

gen e2_carrier = 0
replace e2_carrier = 1 if e2_e2==1
replace e2_carrier = 1 if e2_e3==1
replace e2_carrier = 1 if e2_e4==1

gen e4_carrier = 0
replace e4_carrier = 1 if e2_e4==1
replace e4_carrier = 1 if e3_e4==1
replace e4_carrier = 1 if e4_e4==1


* Inverting e2 homozygous status to have non-e2 homozygotes (e3 or e4 carriers) as risk group:
gen not_e2e2 = 1 if e2_e2!=0
replace not_e2e2 = 0 if e2_e2==1

* Checking there aren't any negative EIDs (withdrawals) in the analytical sample by this point - legacy code from using data locally before UKB-RAP analysis:
count if eid<0

// drop if eid<0


* Defining outcomes:

gen dem = 0
replace dem = 1 if p42019!="" // variables string, rather than numeric here
gen alz = 0
replace alz = 1 if p42021!=""

drop _m

* Adding top-up of case ascertainment from dementia primary care records:
	merge 1:1 eid using "dem_pc_markout.dta"

keep if _m!=2
drop _m

replace dem = 1 if dem_pc==1


merge 1:1 eid using "alz_pc_markout.dta"
	keep if _m!=2
	
	drop _m

replace alz = 1 if alz_pc==1

* Creating ethnicity variable:
tab p21000_i0
tab p21000_i0, nol


gen ethnicity = . if p21000==""

* White category:
replace ethnicity = 1 if p21000=="White" | p21000=="Irish" | p21000=="British" | p21000=="Any other white background"

* Asian category
replace ethnicity = 2 if p21000=="Asian or Asian British" | p21000=="Bangladeshi" | p21000=="Indian" | p21000=="Pakistani" | p21000=="Any other Asian background"

* Black category
replace ethnicity = 3 if p21000=="Black or Black British" | p21000=="Caribbean" | p21000=="African" | p21000=="Any other Black background"

* Mixed category
replace ethnicity = 4 if p21000=="White and Black African" | p21000=="White and Black Caribbean" | p21000=="White and Asian" | p21000=="Mixed" | p21000=="Any other mixed background"

	* Chinese
replace ethnic = 5 if p21000=="Chinese"

	* Other
replace ethnic = 6 if p21000=="Other ethnic group"


label define ethnic_groups 1 "White" 2 "Asian or Asian British" 3 "Black or Black British" 4 "Mixed" 5 "Chinese" 6 "Other"
label values ethnicity ethnic_groups



* Creating binary ethnicity variable so that individuals of rarer ethnicities are not lost from subgroup analyses:
gen bin_eth = 1 if ethnicity==1
replace bin_eth=2 if ethnicity>1 & ethnicity!=.


* Applying standard GWAS exclusions:

gen sex_mismatch =1 if (p22001=="Female" & p31=="Male") | (p22001=="Male" & p31=="Female")

	* aneuploidy
gen aneuploidy = 1 if p22019=="Yes"

	* heterozygosity

tab p22027
gen het_outliers = p22027 

tab sex_mis
tab aneup
tab het_outliers // no obs


* Importing simple relatedness exclusion (randomly pick one from related pairs)
merge 1:1 eid using "rel_exclusions_simple.dta"


keep if _m!=2
drop _m


tab p22021
gen no_kinship = 1 if p22021=="Participant excluded from kinship inference process"

desc sex_mismatch aneuploidy no_kinship simple_rel_exclusions
sum sex_mismatch aneuploidy no_kinship simple_rel_exclusions

egen exclude = rowmax(sex_mismatch aneuploidy no_kinship simple_rel_exclusions) 

drop if exclude == 1


* Creating variable for array type used for genotyping:
tab p22000
gen array = 0 if strpos(p22000, "UKBiLEVE")>0
replace array = 1 if strpos(p22000, "Batch")>0


desc p40000*
destring p40000*, replace

* Creating new sex variable

gen sex = 0 if p31=="Female"
replace sex = 1 if p31=="Male"



* Creating mark-out of analytical sample

mark complete
markout complete not_e2e2 age_baseline sex ethnicity p22009_a1 array
drop if complete!=1


* Identifying follow-up time:
gen ts_40000_i1 = date(p40000_i0, "YMD") 


sum ts_40000_i1 // take largest date value to display below
display %td // shows last record to insert in next line
gen censor_date = date("'insert censor date here'", "DMY") // Current overall censorship date, but reflects incomplete death date coverage (UKB's comprehensive censor date is currently 31st May 2024) and later than censoring of HES records and GP records
* see here for notes about censoring dates: https://biobank.ctsu.ox.ac.uk/showcase/exinfo.cgi?src=Data_providers_and_dates


gen ts_53 = date(p53, "YMD") 


* Maximum follow-up in (death) records:
gen followup = (censor_date - ts_53)/365.25 

sum followup if complete

drop if complete!=1


* Outputting descriptive table:
dtable, cont(age_baseline) fact(sex ethnicity apoe_genotype e2_car e3_car e4_car alz dem) export(ukb_table1_col.xlsx)


* Creating output tables in Excel

putexcel set apoe_addem_RRs_v10, modify 

putexcel A1 = "Exposure (ref group: e2/e2)"
putexcel B1 = "N total"

putexcel C1:H1 = ("AD"), merge hcenter
putexcel I1:N1 = ("All-cause dementia"), merge hcenter

putexcel C2 = "N cases"
putexcel D2 = "Case fraction"
putexcel E2 = "RR"
putexcel F2 = "LCI for RR"
putexcel G2 = "UCI for RR"
putexcel H2 = "Formatted RR"

putexcel I2 = "N cases"
putexcel J2 = "Case fraction"
putexcel K2 = "RR"
putexcel L2 = "LCI for RR"
putexcel M2 = "UCI for RR"
putexcel N2 = "Formatted RR"

putexcel A3 = "e2_e2"


tab apoe_genotype


putexcel A3 = "e2_e2"
putexcel B3 = 1047 //  count from 'tab apoe if complete'


tab alz not_e2e2 if complete
putexcel C3 = 5 // count from 'tab alz not_e2e2'
putexcel C4 = 3376 // count from 'tab alz not_e2e2'


tab dem not_e2e2 if complete
putexcel I3 = 25 // count from 'tab dem not_e2e2' 
putexcel I4 = 7405 // count from 'tab dem not_e2e2'


loc row = 4

foreach var in not_e2e2 {
glm alz `var' age_baseline sex bin_eth p22009_a1 p22009_a2 p22009_a3 p22009_a4 p22009_a5 p22009_a6 p22009_a7 p22009_a8 p22009_a9 p22009_a10 array, family(binomial) link(log) eform
putexcel A`row' = "`var'"
matrix res = r(table)
local N = e(N)

putexcel B`row' = `N'
putexcel E`row' = res[1,1]
putexcel F`row' = res[5,1]
putexcel G`row' = res[6,1]
local ci = "("+string(res[5,1], "%9.2f")+", "+string(res[6,1], "%9.2f")+")"
putexcel H`row' = "`ci'"


glm dem `var' age_baseline sex bin_eth p22009_a1 p22009_a2 p22009_a3 p22009_a4 p22009_a5 p22009_a6 p22009_a7 p22009_a8 p22009_a9 p22009_a10 array, family(binomial) link(log) eform
matrix res = r(table)

putexcel K`row' = res[1,1]
putexcel L`row' = res[5,1]
putexcel M`row' = res[6,1]
local ci = "("+string(res[5,1], "%9.2f")+", "+string(res[6,1], "%9.2f")+")"
putexcel N`row' = "`ci'"
}


glm alz i.apoe sex age_baseline  bin_eth p22009_a1 p22009_a2 p22009_a3 p22009_a4 p22009_a5 p22009_a6 p22009_a7 p22009_a8 p22009_a9 p22009_a10 array, family(binomial) link(log) eform 


matrix res = r(table)
putexcel A6 = "e2_e2"
putexcel A7 = "e2_e3"
putexcel A8 = "e3_e3"
putexcel A9 = "e2_e4"
putexcel A10 = "e3_e4"
putexcel A11 = "e4_e4"


tab apoe_genotype
putexcel B6 = 1047
putexcel B7 = 21018
putexcel B8 = 100892
putexcel B9 = 4208
putexcel B10 = 39974
putexcel B11 = 4032


tab apoe_genotype alz
putexcel C6 = 5
putexcel C7 = 170
putexcel C8 = 1061
putexcel C9 = 81
putexcel C10 = 1608
putexcel C11 = 456


* e2/e3
putexcel E7 = res[1,2]
putexcel F7 = res[5,2]
putexcel G7 = res[6,2]
local ci = "("+string(res[5,2], "%9.2f")+", "+string(res[6,2], "%9.2f")+")"
putexcel H7 = "`ci'"


	* e3_e3
putexcel E8 = res[1,3]
putexcel F8 = res[5,3]
putexcel G8 = res[6,3]
local ci = "("+string(res[5,3], "%9.2f")+", "+string(res[6,3], "%9.2f")+")"
putexcel H8 = "`ci'"

	*  e2/e4
putexcel E9 = res[1,4]
putexcel F9 = res[5,4]
putexcel G9 = res[6,4]
local ci = "("+string(res[5,4], "%9.2f")+", "+string(res[6,4], "%9.2f")+")"
putexcel H9 = "`ci'"

	* e3_e4
putexcel E10 = res[1,5]
putexcel F10 = res[5,5]
putexcel G10 = res[6,5]
local ci = "("+string(res[5,5], "%9.2f")+", "+string(res[6,5], "%9.2f")+")"
putexcel H10 = "`ci'"

	* e4_e4
putexcel E11 = res[1,6]
putexcel F11 = res[5,6]
putexcel G11 = res[6,6]
local ci = "("+string(res[5,6], "%9.2f")+", "+string(res[6,6], "%9.2f")+")"
putexcel H11 = "`ci'"


glm dem i.apoe sex age_baseline  bin_eth p22009_a1 p22009_a2 p22009_a3 p22009_a4 p22009_a5 p22009_a6 p22009_a7 p22009_a8 p22009_a9 p22009_a10 array, family(binomial) link(log) eform 
matrix res = r(table)


tab apoe_genotype dem 
putexcel I6 = 25
putexcel I7 = 508
putexcel I8 = 2893
putexcel I9 = 214
putexcel I10 = 3026
putexcel I11 = 764


* e2/e3
putexcel K7 = res[1,2]
putexcel L7 = res[5,2]
putexcel M7 = res[6,2]
local ci = "("+string(res[5,2], "%9.2f")+", "+string(res[6,2], "%9.2f")+")"
putexcel N7 = "`ci'"


* e3_e3
putexcel K8 = res[1,3]
putexcel L8 = res[5,3]
putexcel M8 = res[6,3]
local ci = "("+string(res[5,3], "%9.2f")+", "+string(res[6,3], "%9.2f")+")"
putexcel N8 = "`ci'"

* e2_e4
putexcel K9 = res[1,4]
putexcel L9 = res[5,4]
putexcel M9 = res[6,4]
local ci = "("+string(res[5,4], "%9.2f")+", "+string(res[6,4], "%9.2f")+")"
putexcel N9 = "`ci'"

* e3_e4 
putexcel K10 = res[1,5]
putexcel L10 = res[5,5]
putexcel M10 = res[6,5]
local ci = "("+string(res[5,5], "%9.2f")+", "+string(res[6,5], "%9.2f")+")"
putexcel N10 = "`ci'"

* e4_e4
putexcel K11 = res[1,6]
putexcel L11 = res[5,6]
putexcel M11 = res[6,6]
local ci = "("+string(res[5,6], "%9.2f")+", "+string(res[6,6], "%9.2f")+")"
putexcel N11 = "`ci'"






