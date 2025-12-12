

* Loading A4 data in:

cd "project directory"
use "2024_07_A4_prerandom_toanalyse.dta", clear

	* NB: see R script from cohort to prepare the pre-randomisation dataset

rename BID bid

	* adding race
	
merge 1:1 bid using "A4_27Aug2024_race_only.dta"
keep if _m==3
drop _m

	* Coding of race variable from data dictionary: 
* PTRACE	PTDEMOG	Participant Demographics	Racial Categories	65	T	1=American Indian or Alaskan Native; 2=Asian; 3=Native Hawaiian or Other Pacific Islander; 4=Black or African American; 5=White; 6=Unknown or Not Reported
	* NB: some entries contain more than one value, presumably admixed individuals reporting more than one group
	
tab ptrace
count if inlist(ptrace, "1:4", "1:4:5", "1:5", "2:5", "4:5") // 27 Mixed
count if inlist(ptrace, "1", "3", "6") // 37 Other

* Recoding  ethnicity variable based on 'race' categories
rename ethnicity hispanic 
replace hispanic = . if hispanic==3

gen ethnicity_binary = 1
replace ethnicity_binary = 0 if ptrace=="5"
replace ethnicity_binary = . if ptrace=="6"

gen subgroups_ethnic = 3 // other group
replace subgroups_ethnic = 0 if ptrace=="5" // white 
replace subgroups_ethnic = 1 if ptrace=="2" // Asian
replace subgroups_ethnic = 2 if ptrace=="4" // Black
replace subgroups_ethnic = . if ptrace=="6" // Other

tab hispanic
tab ethnicity_binary
tab subgroups
tab ethnicity hispanic


	* adding the original eligibility SUVr values used to derive amyloid status:
merge 1:1 bid using "4_PETVADATA_PRV2_05Aug2024.dta"
drop if _m!=3
drop _m

	* NB: variable pmodsuvr is the original SUVr scale used to define amyloid status, whereas the pet_suvr variable is the updated variable
	* In the A4 JAMA Neurology paper, pmodsuvr has been used alongside clinical appraisals to define elevated amyloid status BUT the updated pet_suvr values are listed in table 1, hence there are values of AB+ with SUVr values as low as 0.97



* Coding APOE from string var:

gen e2_e2=0 if apoe!=""
replace e2_e2 = 1 if apoe=="E2/E2"

gen e2_e3=0 if apoe!=""
replace e2_e3 = 1 if apoe=="E2/E3"

gen e2_e4=0 if apoe!=""
replace e2_e4 = 1 if apoe=="E2/E4"


gen e3_e3=0 if apoe!=""
replace e3_e3 = 1 if apoe=="E3/E3"

gen e3_e4=0 if apoe!=""
replace e3_e4 = 1 if apoe=="E3/E4"

gen e4_e4=0 if apoe!=""
replace e4_e4 = 1 if apoe=="E4/E4"


foreach var in e2_e4 e2_e2 e3_e3 e3_e4 e4_e4 e2_e3 {
tab `var'
}


gen apoe_genotype = 1 if e2_e2==1
replace apoe_genotype = 2 if e2_e3==1
replace apoe_genotype = 3 if e2_e4==1
replace apoe_genotype = 4 if e3_e3==1
replace apoe_genotype = 5 if e3_e4==1
replace apoe_genotype = 6 if e4_e4==1

label define apoe 1  "e2_e2" 2 "e2_e3" 3 "e2_e4" 4 "e3_e3" 5 "e3_e4" 6 "e4_e4"
label values apoe_genotype  apoe

tab apoe_genotype

drop apoe

gen e4_carrier = 0 if apoe_genotype!=.
replace e4_carrier = 1 if e2_e4==1 | e3_e4==1 | e4_e4==1

gen e2_carrier = 0 if apoe_genotype!=.
replace e2_carrier = 1 if e2_e4==1 | e2_e3==1 | e2_e2==1

gen e3_carrier=0
replace e3_carrier= 1 if inlist(apoe, 2,4,5)

* Inverting e2 homozygous status to have non-e2 homozygotes (e3 or e4 carriers) as risk group:
gen not_e2e2 = 1 if e2_e2==0
replace not_e2e2 = 0 if e2_e2==1


* Original amyloid status derived to define eligibility (defined by A4 team,  https://jamanetwork.com/journals/jamaneurology/article-abstract/2763540, subsequent advice in emails from Reisa Sperling and Mike Donohue)

gen amypos = 0 if amyloid_status=="negative"
replace amypos = 1 if amyloid_status=="positive"

	* Based on advice from Reisa Sperling, derivation of additional amyloid + status based on 1.15 cutoff using updated SUVr methodology:
	
gen v2_amypos = 0
replace v2_amypos = 1 if pet_suvr>=1.15

* Visualising SUVR data by amyloid status & APOE genotypes:

bysort amypos: sum pmodsuvr // a large number of SUVr values over 1.15 in the negative group due to qualitative reads: 


	* first SUVr definition:
hist pmodsuvr, lcolor(blue%25) fc(blue%25)

twoway (hist pmodsuvr if amypos==0, lcolor(blue%25) fc(blue%25)) ///
(hist pmodsuvr if amypos==1, lcolor(red%25) fc(red%25)), xtitle("Amyloid SUVr") ytitle("Density") legend(order(1 "AB-" 2 "AB+"))

sum pmodsuvr if apoe==1, det
tab amypos if apoe==1
sum pmodsuvr if apoe==1 & amypos==1 
sum pmodsuvr if apoe==4 & amypos==1 
sum pmodsuvr if apoe==6 & amypos==1 


	* Showing mean SUVr values among AB+ group for each set of homozygotes:
twoway (hist pmodsuvr if amypos==0, lcolor(blue%25) fc(blue%25)) ///
(hist pmodsuvr if amypos==1, lcolor(red%25) fc(red%25)) ///
(function 1.19, range(0 6) horizontal  lcolor(green) lpattern(dash) lwidth(thin)) ///
(function 1.29, range(0 6) horizontal  lcolor(orange) lpattern(dash) lwidth(thin)) ///
(function 1.35, range(0 6) horizontal  lcolor(red) lpattern(dash) lwidth(thin)), xtitle("Amyloid SUVr") ytitle("Density") legend(order(1 "AB-" 2 "AB+")) 

	* Showing individual SUVr values for each e2/e2 individual designated as AB+ by original definition
twoway (hist pmodsuvr if amypos==0, lcolor(blue%25) fc(blue%25)) ///
(hist pmodsuvr if amypos==1, lcolor(red%25) fc(red%25) ) ///
(function 1.17, range(0 6) horizontal  lcolor(green) lpattern(dash) lwidth(thin)) ///
(function 1.21, range(0 6) horizontal  lcolor(green) lpattern(dash) lwidth(thin)) 

	* Second SUVr defintion
	
	* overall distribution
twoway (hist pet_suvr, lcolor(red%25) fc(red%25)) ///
(function 1.15, range(0 4) horizontal  lcolor(gray) lpattern(dash) lwidth(thin)) 

sum pet_suvr if e2_e2==1, det 
sum pet_suvr if e3_e3==1 & v2_amypos==1, det 
sum pet_suvr if e4_e4==1 & v2_amypos==1, det 

	* showing e2/e2 positive value on chart
twoway (hist pet_suvr, lcolor(red%25) fc(red%25)) ///
(function 1.15, range(0 4) horizontal  lcolor(black) lpattern(shortdash_dot) lwidth(thin)) ///
(function 1.21, range(0 4) horizontal  lcolor(green) lpattern(dash) lwidth(thin)) 

	* Analytical sample
mark complete 
markout complete v2_amypos apoe_genotype age sex subgroups	// hispanic 


* showing mean SUVrs for each group of homozygotes:

mean pet_suvr if apoe==1 & v2_amypos==1 & complete 
mean pet_suvr if apoe==4 & v2_amypos==1 & complete 
mean pet_suvr if apoe==6 & v2_amypos==1 & complete 


twoway (hist pet_suvr if complete, lcolor(blue%25) fc(blue%25)) ///
(function 1.15, range(0 4) horizontal  lcolor(black) lpattern(shortdash_dot) lwidth(thin)) ///
(function 1.21, range(0 4) horizontal  lcolor(green) lpattern(dash) lwidth(thin)) ///
(function 1.33, range(0 4) horizontal  lcolor(orange) lpattern(dash) lwidth(thin)) ///
(function 1.38, range(0 4) horizontal  lcolor(red) lpattern(dash) lwidth(thin)), legend(order (1 "Density" 2 "Aβ+ threshold" 3 "mean in ε2/ε2 Aβ+" 4 "mean in ε3/ε3 Aβ+" 5 "mean in ε4/ε4 Aβ+")) xtitle("Amyloid SUVr")

graph export "A4_suvr_v4.png", as(png) replace

gen e3e3_v_e4e4 = 1 if apoe==4
replace e3e3_v_e4e4 = 2 if apoe==6

ttest pet_suvr if v2_amypos==1 & complete, by(e3e3_v_e4e4) 

	* For descriptive table:
sum age if complete
tab sex if complete
tab ptrace if complete

tab subgroups if complete, missing

tab apoe if complete

tab e2_carrier if complete
tab e3_carrier if complete
tab e4_carrier if complete

tab v2_amypos if complete
* PTRACE	PTDEMOG	Participant Demographics	Racial Categories	65	T	1=American Indian or Alaskan Native; 2=Asian; 3=Native Hawaiian or Other Pacific Islander; 4=Black or African American; 5=White; 6=Unknown or Not Reported
	* NB: some entries contain more than one value, presumably admixed individuals reporting more than one group?
	
tab ptrace
count if complete & inlist(ptrace, "1:4", "1:4:5", "1:5", "2:5", "4:5") // 27 Mixed
count if complete & inlist(ptrace, "1", "3", "6") // 11 Other - consists of only 1 and 3 in complete-case data, 6s omitted 

tab ptrace if complete	
	
* Automatic output of ORs, RR and PAFs 

* Providing baseline probability of amyloid for RR conversion using formula: OR / (1 - pr + (pr*OR))
	* Formula here:  https://www.bmj.com/content/348/bmj.f7450.short
	bysort apoe: tab amypos
	* 0.08 is baseline probability of AB+ with primary amypos definition in A4 for e2/e2 
	tab not // 0.9944 = fraction exposed to enter in PAF calcs


putexcel set apoe_amyloid_A4_v4, modify sheet("1st AB+ definition") 

putexcel A1 = "Exposure (ref group: e2/e2)"
putexcel B1 = "N total"
putexcel C1 = "N cases"
putexcel D1 = "Case fraction"
putexcel E1:J1 = ("Outcome: amyloid positive, definition 1"), merge hcenter

putexcel E2 = "OR"
putexcel F2 = "95% CI"

putexcel G2 = "RR"
putexcel H2 = "LCI for RR"
putexcel I2 = "UCI for RR"
putexcel J2 = "Formatted RR"


putexcel A3 = "e2_e2"
putexcel B3 = 25
putexcel C3 = 2 // case N based on column for 1 in 'tab apoe amypos if complete'
putexcel C4 = 1300 // case N based on column for 1 in 'tab not_e2e2 amypos if complete'


loc row = 4

foreach var in not_e2e2 {
logistic amypos `var' age sex i.subgroups_ethnic if complete	

putexcel A`row' = "`var'"
matrix res = r(table)
local N = e(N)
local beta = res[1,1]
local lci = res[5,1]
local uci = res[6,1]
local rr = res[1,1] / (1 - 0.08 + (0.08*res[1,1])) // NB: 0.08 is baseline probability of AB+ with original amypos definition in A4 for e2/e2
local rr_lci = res[5,1] / (1 - 0.08 + (0.08*res[5,1])) 
local rr_uci = res[6,1] / (1 - 0.08 + (0.08*res[6,1])) 
 

putexcel B`row' = `N'-25
putexcel E`row' = res[1,1]
local ci = "("+string(res[5,1], "%9.2f")+", "+string(res[6,1], "%9.2f")+")"
putexcel F`row' = "`ci'"
putexcel G`row' = `rr'
putexcel H`row' = `rr_lci'
putexcel I`row' = "`rr_uci'"
local rr_ci = "("+string(`rr_lci', "%9.2f")+", "+string(`rr_uci', "%9.2f")+")"
putexcel J`row' = "`rr_ci'"
}

	* Outputting results from individual genotypes for amyloid
logistic amypos i.apoe_genotype age sex i.subgroups_ethnic if complete


matrix res = r(table)
putexcel A6 = "e2_e2"
putexcel A7 = "e2_e3"
putexcel A8 = "e2_e4"
putexcel A9 = "e3_e3"
putexcel A10 = "e3_e4"
putexcel A11 = "e4_e4"

	// NB: Ns below generated by 'tab apoe_genotype if complete'
putexcel B6 = 25
putexcel B7 = 446
putexcel B8 = 115
putexcel B9 = 2400
putexcel B10 = 1291
putexcel B11 = 138

putexcel C6 = 2 // following case Ns based on column for 1 in 'tab apoe amypos if complete'
putexcel C7 = 67
putexcel C8 = 42
putexcel C9 = 477
putexcel C10 = 609
putexcel C11 = 105


loc row = 7
loc pos = 2

	* e2/e3
putexcel E`row' = res[1,`pos']
local ci = "("+string(res[5,`pos'], "%9.2f")+", "+string(res[6,`pos'], "%9.2f")+")"
putexcel F`row' = "`ci'"

local rr = res[1,`pos'] / (1 - 0.08 + (0.08*res[1,`pos'])) 
local rr_lci = res[5,`pos'] / (1 - 0.08 + (0.08*res[5,`pos'])) 
local rr_uci = res[6,`pos'] / (1 - 0.08 + (0.08*res[6,`pos'])) 
 
putexcel G`row' = `rr'
putexcel H`row' = `rr_lci'
putexcel I`row' = `rr_uci'
local rr_ci = "("+string(`rr_lci', "%9.2f")+", "+string(`rr_uci', "%9.2f")+")"
putexcel J`row' = "`rr_ci'"

loc row = `row' + 1
loc pos = `pos' + 1
 
	* e2/e4
putexcel E`row' = res[1,`pos']
local ci = "("+string(res[5,`pos'], "%9.2f")+", "+string(res[6,`pos'], "%9.2f")+")"
putexcel F`row' = "`ci'"

local rr = res[1,`pos'] / (1 - 0.08 + (0.08*res[1,`pos'])) 
local rr_lci = res[5,`pos'] / (1 - 0.08 + (0.08*res[5,`pos'])) 
local rr_uci = res[6,`pos'] / (1 - 0.08 + (0.08*res[6,`pos'])) 
 
putexcel G`row' = `rr'
putexcel H`row' = `rr_lci'
putexcel I`row' = `rr_uci'
local rr_ci = "("+string(`rr_lci', "%9.2f")+", "+string(`rr_uci', "%9.2f")+")"
putexcel J`row' = "`rr_ci'"

loc row = `row' + 1
loc pos = `pos' + 1

 
	* e3_e3
putexcel E`row' = res[1,`pos']
local ci = "("+string(res[5,`pos'], "%9.2f")+", "+string(res[6,`pos'], "%9.2f")+")"
putexcel F`row' = "`ci'"

local rr = res[1,`pos'] / (1 - 0.08 + (0.08*res[1,`pos'])) 
local rr_lci = res[5,`pos'] / (1 - 0.08 + (0.08*res[5,`pos'])) 
local rr_uci = res[6,`pos'] / (1 - 0.08 + (0.08*res[6,`pos'])) 
 
putexcel G`row' = `rr'
putexcel H`row' = `rr_lci'
putexcel I`row' = `rr_uci'
local rr_ci = "("+string(`rr_lci', "%9.2f")+", "+string(`rr_uci', "%9.2f")+")"
putexcel J`row' = "`rr_ci'"

loc row = `row' + 1
loc pos = `pos' + 1

 
	* e3_e4
putexcel E`row' = res[1,`pos']
local ci = "("+string(res[5,`pos'], "%9.2f")+", "+string(res[6,`pos'], "%9.2f")+")"
putexcel F`row' = "`ci'"

local rr = res[1,`pos'] / (1 - 0.08 + (0.08*res[1,`pos'])) 
local rr_lci = res[5,`pos'] / (1 - 0.08 + (0.08*res[5,`pos'])) 
local rr_uci = res[6,`pos'] / (1 - 0.08 + (0.08*res[6,`pos'])) 
 
putexcel G`row' = `rr'
putexcel H`row' = `rr_lci'
putexcel I`row' = `rr_uci'
local rr_ci = "("+string(`rr_lci', "%9.2f")+", "+string(`rr_uci', "%9.2f")+")"
putexcel J`row' = "`rr_ci'"

loc row = `row' + 1
loc pos = `pos' + 1

	* e4_e4
putexcel E`row' = res[1,`pos']
local ci = "("+string(res[5,`pos'], "%9.2f")+", "+string(res[6,`pos'], "%9.2f")+")"
putexcel F`row' = "`ci'"

local rr = res[1,`pos'] / (1 - 0.08 + (0.08*res[1,`pos'])) 
local rr_lci = res[5,`pos'] / (1 - 0.08 + (0.08*res[5,`pos'])) 
local rr_uci = res[6,`pos'] / (1 - 0.08 + (0.08*res[6,`pos'])) 
 
putexcel G`row' = `rr'
putexcel H`row' = `rr_lci'
putexcel I`row' = `rr_uci'
local rr_ci = "("+string(`rr_lci', "%9.2f")+", "+string(`rr_uci', "%9.2f")+")"
putexcel J`row' = "`rr_ci'"





*** Repeating results generation for outcome definition #2


putexcel set apoe_amyloid_A4_v4, modify sheet("2nd AB+ definition") 

putexcel A1 = "Exposure (ref group: e2/e2)"
putexcel B1 = "N total"
putexcel C1 = "N cases"
putexcel D1 = "Case fraction"
putexcel E1:J1 = ("Outcome: amyloid positive, definition 2"), merge hcenter

putexcel E2 = "OR"
putexcel F2 = "95% CI"

putexcel G2 = "RR"
putexcel H2 = "LCI for RR"
putexcel I2 = "UCI for RR"
putexcel J2 = "Formatted RR"


putexcel A3 = "e2_e2"
putexcel B3 = 25
putexcel C3 = 1 // case N based on column for 1 in 'tab apoe v2_amypos if complete'
putexcel C4 = 1202 // case N based on column for 1 in 'tab not_e2e2 v2_amypos if complete'


loc row = 4

foreach var in not_e2e2 {
logistic v2_amypos `var' age sex i.subgroups_ethnic if complete

putexcel A`row' = "`var'"
matrix res = r(table)
local N = e(N)
local beta = res[1,1]
local lci = res[5,1]
local uci = res[6,1]
local rr = res[1,1] / (1 - 0.04 + (0.04*res[1,1])) // NB: 0.04 is baseline probability of AB+ with updated analytical amypos definition in A4 for e2/e2
local rr_lci = res[5,1] / (1 - 0.04 + (0.04*res[5,1])) 
local rr_uci = res[6,1] / (1 - 0.04 + (0.04*res[6,1])) 
 

putexcel B`row' = `N'-25
putexcel E`row' = res[1,1]
local ci = "("+string(res[5,1], "%9.2f")+", "+string(res[6,1], "%9.2f")+")"
putexcel F`row' = "`ci'"
putexcel G`row' = `rr'
putexcel H`row' = `rr_lci'
putexcel I`row' = "`rr_uci'"
local rr_ci = "("+string(`rr_lci', "%9.2f")+", "+string(`rr_uci', "%9.2f")+")"
putexcel J`row' = "`rr_ci'"
}

	* Outputting results from individual genotypes for amyloid
logistic v2_amypos i.apoe_genotype age sex i.subgroups_ethnic if complete


matrix res = r(table)
putexcel A6 = "e2_e2"
putexcel A7 = "e2_e3"
putexcel A8 = "e2_e4"
putexcel A9 = "e3_e3"
putexcel A10 = "e3_e4"
putexcel A11 = "e4_e4"

	// NB: Ns below generated by 'tab apoe_genotype if complete'
putexcel B6 = 25
putexcel B7 = 446
putexcel B8 = 115
putexcel B9 = 2400
putexcel B10 = 1291
putexcel B11 = 138

putexcel C6 = 1 // following case Ns based on column for 1 in 'tab apoe v2_amypos if complete'
putexcel C7 = 60
putexcel C8 = 37
putexcel C9 = 409
putexcel C10 = 591
putexcel C11 = 105


loc row = 7
loc pos = 2

	* e2/e3
putexcel E`row' = res[1,`pos']
local ci = "("+string(res[5,`pos'], "%9.2f")+", "+string(res[6,`pos'], "%9.2f")+")"
putexcel F`row' = "`ci'"

local rr = res[1,`pos'] / (1 - 0.04 + (0.04*res[1,`pos'])) 
local rr_lci = res[5,`pos'] / (1 - 0.04 + (0.04*res[5,`pos'])) 
local rr_uci = res[6,`pos'] / (1 - 0.04 + (0.04*res[6,`pos'])) 
 
putexcel G`row' = `rr'
putexcel H`row' = `rr_lci'
putexcel I`row' = `rr_uci'
local rr_ci = "("+string(`rr_lci', "%9.2f")+", "+string(`rr_uci', "%9.2f")+")"
putexcel J`row' = "`rr_ci'"

loc row = `row' + 1
loc pos = `pos' + 1
 
	* e2/e4
putexcel E`row' = res[1,`pos']
local ci = "("+string(res[5,`pos'], "%9.2f")+", "+string(res[6,`pos'], "%9.2f")+")"
putexcel F`row' = "`ci'"

local rr = res[1,`pos'] / (1 - 0.04 + (0.04*res[1,`pos'])) 
local rr_lci = res[5,`pos'] / (1 - 0.04 + (0.04*res[5,`pos'])) 
local rr_uci = res[6,`pos'] / (1 - 0.04 + (0.04*res[6,`pos'])) 
 
putexcel G`row' = `rr'
putexcel H`row' = `rr_lci'
putexcel I`row' = `rr_uci'
local rr_ci = "("+string(`rr_lci', "%9.2f")+", "+string(`rr_uci', "%9.2f")+")"
putexcel J`row' = "`rr_ci'"

loc row = `row' + 1
loc pos = `pos' + 1

 
	* e3_e3
putexcel E`row' = res[1,`pos']
local ci = "("+string(res[5,`pos'], "%9.2f")+", "+string(res[6,`pos'], "%9.2f")+")"
putexcel F`row' = "`ci'"

local rr = res[1,`pos'] / (1 - 0.04 + (0.04*res[1,`pos'])) 
local rr_lci = res[5,`pos'] / (1 - 0.04 + (0.04*res[5,`pos'])) 
local rr_uci = res[6,`pos'] / (1 - 0.04 + (0.04*res[6,`pos'])) 
 
putexcel G`row' = `rr'
putexcel H`row' = `rr_lci'
putexcel I`row' = `rr_uci'
local rr_ci = "("+string(`rr_lci', "%9.2f")+", "+string(`rr_uci', "%9.2f")+")"
putexcel J`row' = "`rr_ci'"

loc row = `row' + 1
loc pos = `pos' + 1

 
	* e3_e4
putexcel E`row' = res[1,`pos']
local ci = "("+string(res[5,`pos'], "%9.2f")+", "+string(res[6,`pos'], "%9.2f")+")"
putexcel F`row' = "`ci'"

local rr = res[1,`pos'] / (1 - 0.04 + (0.04*res[1,`pos'])) 
local rr_lci = res[5,`pos'] / (1 - 0.04 + (0.04*res[5,`pos'])) 
local rr_uci = res[6,`pos'] / (1 - 0.04 + (0.04*res[6,`pos'])) 
 
putexcel G`row' = `rr'
putexcel H`row' = `rr_lci'
putexcel I`row' = `rr_uci'
local rr_ci = "("+string(`rr_lci', "%9.2f")+", "+string(`rr_uci', "%9.2f")+")"
putexcel J`row' = "`rr_ci'"

loc row = `row' + 1
loc pos = `pos' + 1

	* e4_e4
putexcel E`row' = res[1,`pos']
local ci = "("+string(res[5,`pos'], "%9.2f")+", "+string(res[6,`pos'], "%9.2f")+")"
putexcel F`row' = "`ci'"

local rr = res[1,`pos'] / (1 - 0.04 + (0.04*res[1,`pos'])) 
local rr_lci = res[5,`pos'] / (1 - 0.04 + (0.04*res[5,`pos'])) 
local rr_uci = res[6,`pos'] / (1 - 0.04 + (0.04*res[6,`pos'])) 
 
putexcel G`row' = `rr'
putexcel H`row' = `rr_lci'
putexcel I`row' = `rr_uci'
local rr_ci = "("+string(`rr_lci', "%9.2f")+", "+string(`rr_uci', "%9.2f")+")"
putexcel J`row' = "`rr_ci'"
