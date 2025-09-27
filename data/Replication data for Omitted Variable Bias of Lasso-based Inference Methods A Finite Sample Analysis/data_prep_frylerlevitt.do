*****************************************************************************************
*** Data preparation application Fryer and Levitt (2013, AER)
*** Authors: Kaspar Wuthrich and Ying Zhu
*** DISCLAIMER: This software is provided "as is" without warranty of any kind, expressed or implied. 
*** Questions/error reports: kwuthrich@ucsd.edu
*****************************************************************************************


clear all
set more off

cd "/Users/kasparwuthrich/Dropbox/research/Kaspar_Ying_research/submission/ReStat/final submission/replication_package_final"

use "cpp.dta", clear

*****************************************************************************************
*****************************************************************************************
*** Data preparation as in Fryer and Levitt
/*
COPYRIGHT 2013 American Economic Association

=================================================================
Modified BSD License
=================================================================

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
*****************************************************************************************
*****************************************************************************************

******************************************
*** Fix Up some variables
******************************************

*** Race
gen race = race_8months
replace race = race_7years if race == .
replace race = race_3years if race == .
replace race = race_4years if race == .

*** Sex
gen sex = sex_8months
replace sex = sex_7years if sex == .
replace sex = sex_3years if sex == .
replace sex = sex_4years if sex == .

*** Female
gen female = (sex == 2)

*** Examiner IDs
gen examiner_id_8months = institution_8months + examiner_8months
gen examiner_id_4years = institution_4years + examiner_4years
gen examiner_id_7years = institution_7years + examiner_7years

*** Race
gen white = (race == 1)
gen black = (race == 2)
gen hispanic = (race == 4)
gen other = (race == 3 | race == 8 | race == 9)

*** Multiple births
egen multiple = max(fetus_number), by(mother_child_id)
gen singleton = (multiple == 0 | multiple == 1)
gen twin = (multiple == 2)
gen high_order_multiple = (multiple == 3)

*** Birth Weight
replace birth_weight = "" if birth_weight == "9999"
gen birth_weight_pds = substr(birth_weight, 1, 2)
gen birth_weight_ounce = substr(birth_weight, 3, 2)
destring(birth_weight_pds), replace
destring(birth_weight_ounce), replace
gen weight_at_birth = birth_weight_pds * 16 + birth_weight_ounce
gen w_less_1500 = (weight_at_birth <= 53)
gen w_1500_2500 = (weight_at_birth > 53 & weight_at_birth <= 88)
gen w_2500_3500 = (weight_at_birth > 88 & weight_at_birth <= 124) 
gen w_3500_more = (weight_at_birth > 124 & weight_at_birth != .) 

*** Schooling
gen dad_hs_dropout = (dad_highest_grade < 12)
gen dad_hs_grad = (dad_highest_grade == 12)
gen dad_some_college = (dad_highest_grade > 12 & dad_highest_grade < 16)
gen dad_college_plus = (dad_highest_grade == 16 | dad_highest_grade == 17 | dad_highest_grade == 18)
*
gen mom_hs_dropout = (mom_highest_grade < 12)
gen mom_hs_grad = (mom_highest_grade == 12)
gen mom_some_college = (mom_highest_grade > 12 & mom_highest_grade < 16)
gen mom_college_plus = (mom_highest_grade == 16 | mom_highest_grade == 17 | mom_highest_grade == 18)

*** Occupation
gen dad_no_occupation = (dad_occupation == 0 | dad_occupation == 1)
gen dad_professional = (dad_occupation == 10 | dad_occupation == 12 | dad_occupation == 20 | dad_occupation == 82)
gen dad_non_professional = (dad_occupation == 30 | dad_occupation == 40 | dad_occupation == 50 | dad_occupation == 60 | dad_occupation == 70 | dad_occupation == 72 | dad_occupation == 80)
*
gen mom_no_occupation = (mom_occupation == 0 | mom_occupation == 1)
gen mom_professional = (mom_occupation == 10 | mom_occupation == 12 | mom_occupation == 20 | mom_occupation == 82)
gen mom_non_professional = (mom_occupation == 30 | mom_occupation == 40 | mom_occupation == 50 | mom_occupation == 60 | mom_occupation == 70 | mom_occupation == 72 | mom_occupation == 80)

*** Income
gen inc_less_500 = (income == 0 | income == 5 | income == 15)
gen inc_500_1000 = (income == 25 | income == 35)
gen inc_1000_1500 = (income == 45 | income == 55)
gen inc_1500_2000 = (income == 65 | income == 75)
gen inc_2000_2500 = (income == 85 | income == 95)
gen inc_2500_plus = (income == 96)

*** Parents Combo
gen both_bio_parents = (is_real_mom == 1 & (husband_home == 1 | husband_home == 2))
gen other_parent_config = (both_bio_parents == 0)

*** Premature
replace weeks_gestation = . if (weeks_gestation == 88 | weeks_gestation ==  99)
gen weeks_premature = 37 - weeks_gestation
gen weeks_premature_0 = (weeks_premature <= 0)
gen weeks_premature_1 = (weeks_premature == 1) 
gen weeks_premature_2 = (weeks_premature == 2) 
gen weeks_premature_3 = (weeks_premature == 3)
gen weeks_premature_4 = (weeks_premature == 4) 
gen weeks_premature_5 = (weeks_premature == 5) 
gen weeks_premature_6 = (weeks_premature == 6) 
gen weeks_premature_7 = (weeks_premature == 7) 
gen weeks_premature_8 = (weeks_premature == 8) 
gen weeks_premature_9 = (weeks_premature == 9) 
gen weeks_premature_10 = (weeks_premature == 10) 
gen weeks_premature_11 = (weeks_premature >= 11 & weeks_premature != .) 

*** Parental Score
gen mother_negative = (express_affection == 1 | evaluation_child == 1 | handling_child == 1 | managment_child == 1 | reaction_child_need == 1 | reaction_child_test == 1 | focus_on_atention == 1 | child_appearance == 1)
gen mother_indifferent = (express_affection == 2 | evaluation_child == 2 | handling_child == 2 | managment_child == 2 | reaction_child_need == 2 | reaction_child_test == 2 | focus_on_atention == 2 | child_appearance == 2)
gen mother_accepting = (express_affection == 3 | evaluation_child == 3 | handling_child == 3 | managment_child == 3 | reaction_child_need == 3 | reaction_child_test == 3 | focus_on_atention == 3 | child_appearance == 3)
gen mother_attentive = (express_affection == 4 | evaluation_child == 4 | handling_child == 4 | managment_child == 4 | reaction_child_need == 4 | reaction_child_test == 4 | focus_on_atention == 4 | child_appearance == 4)
gen mother_over_caring = (express_affection == 5 | evaluation_child == 5 | handling_child == 5 | managment_child == 5 | reaction_child_need == 5 | reaction_child_test == 5 | focus_on_atention == 5 | child_appearance == 5)
gen mother_other = (express_affection == 6 | evaluation_child == 6 | handling_child == 6 | managment_child == 6 | reaction_child_need == 6 | reaction_child_test == 6 | focus_on_atention == 6 | child_appearance == 6)
foreach var of varlist express_affection evaluation_child handling_child managment_child reaction_child_need reaction_child_test focus_on_atention child_appearance{
  replace `var' = . if (`var' == 8 | `var' == 9)
}

**** AGES
gen age_8months_1 = (age_8months == 1)
gen age_8months_2 = (age_8months == 2)
gen age_8months_3 = (age_8months == 3)
gen age_8months_4 = (age_8months == 4)
gen age_8months_5 = (age_8months == 5)
gen miss_age8months = (age_8months == 9)

gen miss_age_4years = (age_4years == 19)
replace age_4years = 0 if age_4years == 19
tostring(age_4years), replace
gen age_4years_years = substr(age_4years, 1, 1)
gen age_4years_months = substr(age_4years, 2, 2)
rename age_4years age_4years_original
destring(age_4years_years), replace
destring(age_4years_months), replace
gen age_4years = age_4years_years * 12 + age_4years_months
replace age_4years = 0 if age_4years == .
replace age_4years = 45 if age_4years < 45 & age_4years > 0 
replace age_4years = 55 if age_4years > 55
drop age_4years_years age_4years_months

replace age_7years = 500 if age_7years == 5
replace age_7years = 600 if age_7years == 6
replace age_7years = 700 if age_7years == 7
replace age_7years = 800 if age_7years == 8
replace age_7years = 900 if age_7years == 9
replace age_7years = . if age_7years < 606 & age_7years > 800
gen miss_age_7years = (age_7years == .)
replace age_7years = 0 if age_7years == .
tostring(age_7years), replace
gen age_7years_years = substr(age_7years, 1, 1)
gen age_7years_months = substr(age_7years, 2, 2)
destring(age_7years_years), replace
destring(age_7years_months), replace
rename age_7years age_7years_original
gen age_7years = age_7years_years * 12 + age_7years_months
replace age_7years = 0 if age_7years == .
replace age_7years = 81 if age_7years < 81 & age_7years > 0
replace age_7years = 99 if age_7years > 99
drop age_7years_years age_7years_months


*** replace scores if missing
replace mental_score_8months = . if mental_score_8months == 999 | age_8months == 9
replace iq_4years = . if iq_4years == 999 
replace fullscale_iq_7years = . if (fullscale_iq_7years == 999 | fullscale_iq_7years == 777 | fullscale_iq_7years == 888) | age_7years == 0
gen miss_raw_score_mother = (raw_score_mother == . | raw_score_mother == 88)
replace raw_score_mother = . if raw_score_mother == 88
* we drop 8-month old score that are less than -10SD
egen stand_mental_score_8months = std(mental_score_8months)
replace mental_score_8months = . if stand_mental_score_8months < -10
drop stand_mental_score_8months

** Drop if all scores are missing
keep if (mental_score_8months != . & iq_4years != . & fullscale_iq_7years != . & examiner_4years != "") 

** Standardize Scores
foreach var of varlist raw_score_mother mental_score_8months iq_4years fullscale_iq_7years{
	egen stand_`var' = std(`var')
}
replace stand_raw_score_mother = 0 if  miss_raw_score_mother == 1

*** Siblings
gen siblings_0 = (siblings == 0)
gen siblings_1 = (siblings == 1)
gen siblings_2 = (siblings == 2)
gen siblings_3 = (siblings == 3)
gen siblings_4 = (siblings == 4)
gen siblings_5 = (siblings == 5)
gen siblings_6_plus = (siblings >= 6 & siblings != .)

*** Replace age_mom
gen miss_age_mom = (age_mom == 99 | age_mom == .)
replace age_mom = 0 if miss_age_mom == 1
gen age_mom_2 = age_mom * age_mom
gen age_mom_3 = age_mom_2 * age_mom
gen age_mom_4 = age_mom_3 * age_mom / 1000
gen age_mom_5 = age_mom_4 * age_mom / 1000

*** Missing Dummies
gen miss_income = (income == . | income == 99)
gen miss_dad_occupation = (dad_occupation == . | dad_occupation == 99)
gen miss_mom_occupation = (mom_occupation == . | mom_occupation == 99) 
gen miss_mom_education = (mom_highest_grade == . | mom_highest_grade == 77 | mom_highest_grade == 99)
gen miss_dad_education = (dad_highest_grade == . | dad_highest_grade == 77 | dad_highest_grade == 99)
gen miss_birthweight = (weight_at_birth == .)
gen miss_multiple = (multiple == .)
gen miss_parents  = (husband_home == . | husband_home == 9) | (is_real_mom == . | is_real_mom == 9)
gen miss_siblings = (siblings == . | siblings == 99)
gen miss_premature = (weeks_gestation == .)
gen miss_parental_score = (express_affection == . & evaluation_child == . & handling_child == . & managment_child == . & reaction_child_need == . & reaction_child_test == . & focus_on_atention == . & child_appearance == .)
gen miss_mom_age = age_mom == 0
gen moms_age = age_mom
replace moms_age = . if age_mom == 0
gen age_4year_olds = age_4years
replace age_4year_olds = . if age_4years == 0
gen age_7year_olds = age_7years
replace age_7year_olds = . if age_7years == 0


*** local groupings
global basic_8months "i.age_8months female"
global basic_4years "i.age_4years female"
global basic_7years "i.age_7years female"
global ses "dad_hs_dropout dad_hs_grad dad_some_college dad_college_plus dad_no_occupation dad_professional dad_non_professional mom_hs_dropout mom_hs_grad mom_some_college mom_college_plus mom_no_occupation mom_professional mom_non_professional inc_less_500 inc_500_1000 inc_1000_1500 inc_1500_2000 inc_2000_2500 inc_2500_plus"
global home_environment "siblings_0 siblings_1 siblings_2 siblings_3 siblings_4 siblings_5 siblings_6_plus  both_bio_parents  age_mom age_mom_2 age_mom_3 age_mom_4 age_mom_5 miss_age_mom mother_indifferent mother_accepting mother_attentive mother_over_caring mother_other miss_parental_score"
global prenatal "w_less_1500 w_1500_2500 w_2500_3500 w_3500_more weeks_premature_*  singleton twin high_order_multiple "
global mom_score "stand_raw_score_mother miss_raw_score_mother"

* omitting other_parent_config miss_premature miss_birthweight miss_multiple miss_siblings miss_income miss_dad_occupation miss_mom_occupation miss_age_7years, age4 <= 45 age missing8

******************************************
*** Summary Statistics 
*****************************************
gen premature = weeks_premature < 0
foreach condition in "" " if white == 1" " if black == 1" " if hispanic == 1" " if other == 1"{
  di "*****************************"
  di "`condition'"
  count `condition'
  summ stand_mental_score_8months stand_iq_4years stand_fullscale_iq_7years `condition'
  summ age_8months_* age_4year_olds age_7year_olds female `condition'
  summ miss_dad_education miss_dad_occupation miss_mom_education miss_mom_occupation miss_income miss_siblings miss_parents miss_mom_age miss_parental_score miss_raw_score_mother miss_birthweight miss_premature miss_multiple  `condition'
}

foreach condition in " " " & white == 1" " & black == 1" " & hispanic == 1" " & other == 1"{
  di "*****************************"
  di "`condition'"
  summ dad_hs_dropout dad_hs_grad dad_some_college dad_college_plus if miss_dad_education != 1 `condition'
  summ dad_no_occupation dad_professional dad_non_professional if miss_dad_occupation != 1 `condition'
  summ mom_hs_dropout mom_hs_grad mom_some_college mom_college_plus if miss_mom_education != 1 `condition'
  summ mom_no_occupation mom_professional mom_non_professional if miss_mom_occupation != 1 `condition'
  summ inc_less_500 inc_500_1000 inc_1000_1500 inc_1500_2000 inc_2000_2500 inc_2500_plus if miss_income != 1 `condition'
  summ siblings if miss_siblings != 1 `condition'
  summ both_bio_parents other_parent_config if miss_parents != 1 `condition'
  summ moms_age if miss_mom_age != 1 `condition'
  summ mother_indifferent mother_accepting mother_attentive mother_over_caring mother_other  if miss_parental_score != 1 `condition'
  summ stand_raw_score_mother if miss_raw_score_mother != 1 `condition'
  summ w_less_1500 w_1500_2500 w_2500_3500 w_3500_more  if miss_birthweight != 1 `condition'
  summ premature if miss_premature != 1 `condition'
  summ weeks_premature if premature > 0 `condition'
  summ singleton twin high_order_multiple if miss_multiple != 1 `condition'
}

*** Number of missing values 
gen missing_sum = miss_dad_education + miss_dad_occupation + miss_mom_education + miss_mom_occupation + miss_income + miss_siblings + miss_parents + miss_mom_age + miss_parental_score + miss_raw_score_mother + miss_birthweight + miss_premature + miss_multiple 
gen missing_total = "0" if missing_sum == 0
replace missing_total = "1-2" if missing_sum == 1 | missing_sum == 2
replace missing_total = "3+" if missing_sum > 2
tab missing_total 
tab missing_total if white == 1
tab missing_total if black == 1
tab missing_total if hispanic == 1
tab missing_total if other == 1


*****************************************************************************************
*****************************************************************************************
*** End of the data preparation in Fryer and Levitt
*****************************************************************************************
*****************************************************************************************

******************************************
***  Replication of results in Table 3
******************************************

*** 7 YEARS
gen interviewer_7years = institution_7years + examiner_7years

xi: reg stand_fullscale_iq_7years black hispanic other
xi: areg stand_fullscale_iq_7years black hispanic other $basic_7years $ses $home_environment $prenatal, absorb(interviewer_7years)

******************************************
***  Prepare dataset for analysis
******************************************

keep if e(sample)

* Focus on black-white gap
drop if hispanic == 1
drop if other ==1

tab age_8months, gen(dage8m_)
tab age_4years, gen(dage4y_)
tab age_7years, gen(dage7y_)

save data_test_scores.dta, replace

* Test regression

reg stand_fullscale_iq_7years black $basic_7years $ses $home_environment $prenatal,robust

* Regression for R

reg stand_fullscale_iq_7years black dage7y_1 dage7y_2 dage7y_3 dage7y_4 dage7y_5 dage7y_6 dage7y_7 dage7y_8 dage7y_9 dage7y_10 dage7y_11 dage7y_12 dage7y_13 dage7y_14 dage7y_15 dage7y_16 dage7y_17 dage7y_18 female dad_hs_dropout dad_hs_grad dad_some_college dad_college_plus dad_no_occupation dad_professional dad_non_professional mom_hs_dropout mom_hs_grad mom_some_college mom_college_plus mom_no_occupation mom_professional mom_non_professional inc_less_500 inc_500_1000 inc_1000_1500 inc_1500_2000 inc_2000_2500 inc_2500_plus siblings_0 siblings_1 siblings_2 siblings_3 siblings_4 siblings_5 siblings_6_plus  both_bio_parents  age_mom age_mom_2 age_mom_3 age_mom_4 age_mom_5 miss_age_mom mother_indifferent mother_accepting mother_attentive mother_over_caring mother_other miss_parental_score w_less_1500 w_1500_2500 w_2500_3500 w_3500_more weeks_premature_0-weeks_premature_11 singleton twin high_order_multiple, robust


reg stand_fullscale_iq_7years black, r


