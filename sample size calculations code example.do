
...

svyset hv001 [weight=new_wgt], singleunit(centered) strata(hv022)

/* GENERATING STRATA VARIABLE */

egen strata_grouped = group(hv024 hv025 hv104), label

***** create dummy for relevant strata groups

	label list strata_grouped

	global relevant_groups 17 18 19 20     65 66 67 68     101 102 103 104     109 110 111 112     129 130 131 132

	capture drop relevant_group
	gen relevant_group = .
	foreach group of numlist $relevant_groups {
		replace relevant_group = 1 if strata_grouped == `group' & relevant_group == .
	}

***** save the total number of relevant groups	
	
	tab strata_grouped if relevant_group == 1
	return list
	global max_groups = `r(r)'
	di $max_groups

/* CALCULATE DESIGN EFFECTS */

***** obtain mean/proportion for each group

	svy , subpop(if age_schl>=6 & age_schl<=10 & relevant_group == 1) : mean enrolled_primary , over(strata_grouped)

	matrix effects=r(table)

	foreach n of numlist 1/ $max_groups {
		global original_mean`n' = effects[1,`n']
	}
	
	di $original_mean20

***** obtain deff for each group

	estat effects , srssubpop

	foreach n of numlist 1/ $max_groups {
		global original_deff`n' = r(deffsub)[1,`n']
	}
	
	di $original_deff20

***** obtain standard deviation for each group

	estat sd , srssubpop

	foreach n of numlist 1/ $max_groups {
		global original_sd`n' = r(sd)[1,`n']
	}
	
	di $original_sd20
	
***** calculate infinite n SRS for each group

	global tt = 1.96 * 1.96 // 95% certainty ... adjustable!
	global precision = .05 // + or - 5% ... adjustable!
	
	foreach n of numlist 1/ $max_groups {
		local variance = ${original_sd`n'} * ${original_sd`n'}
		local E = ( $precision * ${original_mean`n'} )
		global insrs`n' = ( $tt * `variance' ) / (`E' * `E')
	}
	
	di $insrs20
	
***** calculate finite n SRS

	// globals with population for each strata group from https://www.niti.gov.in/niti/content/population-number-male-female-rural-urban
	
	local n = 0
	foreach pop of numlist 6201496 5528113 47983851 44091177 7617584 8314587 8403706 9051800 3616819 3379305 17584859 17366375 5548353 4839083 9086466 8230334 23551760 20918695 81044655 74066367 {
		local n = `n' + 1
		global pop`n' = `pop'
	}
	
	di $pop20

	foreach n of numlist 1/ $max_groups {
		global fnsrs`n' = ( ${insrs`n'} ) / (1 + ( ${insrs`n'} / ${pop`n'} ) )
	}
	
	di $fnsrs20

***** obtain original m for each group

	capture drop n
	gen n = 1 if enrolled_primary !=.
	preserve
		collapse (sum) n, by(strata_grouped hv001)
		local n = 0
		foreach group of numlist $relevant_groups {
			local n = `n' + 1
			sum n if strata_grouped == `group'
			quietly return list
			global original_m`n' = `r(mean)'
		}
	restore
	
	di $original_m20

***** calculate rho using DEFF-1/m-1

	foreach n of numlist 1/ $max_groups {
		global rho`n' = ( ${original_deff`n'} - 1 ) / ( ${original_m`n'} - 1 )
	}
	di $rho20
	
***** calculate DEFF for our survey

	global new_m = 30 // cluster size ... adjustable!
	
	foreach n of numlist 1/ $max_groups {
		global new_deff`n' = 1 + ( ${rho`n'} * ( $new_m - 1 ) )
	}

	di $new_deff20
	
***** apply design effects to finite population

	capture drop strata_grouped_s
	decode strata_grouped , gen(strata_grouped_s)
	local n = 0
	foreach group of numlist $relevant_groups {
		local n = `n' + 1
		levelsof strata_grouped_s if strata_grouped == `group', local(desc) clean
		global label`n' "`desc' (strata `group'/order `n')"
	}

	foreach n of numlist 1/ $max_groups {
		global n`n' =  ${fnsrs`n'} * ${new_deff`n'}
		global nrounded`n' = ceil( ${n`n'})
		di "Sample size for ${label`n'} is ${nrounded`n'}"
	}
	
	di $nrounded20
	
***** calculate total sample size needed if some portion of calls are not picked up
	
	global nonresponse = .4 // 40% of calls are not picked up ... adjustable!
	
	foreach n of numlist 1/ $max_groups {
		global nNR`n' = (  ${n`n'} * $nonresponse ) +  ${n`n'}
		global nNRrounded`n' = ceil( ${nNR`n'})
		di "Sample size for ${label`n'}, accounting for non-response, is ${nNRrounded`n'}"	
	}
	
	di $nNRrounded20
	
***** calculate total sample size needed overall

	global total_sample = 0
	global total_sample_NR = 0
	foreach n of numlist 1/ $max_groups {
		global total_sample = $total_sample + ${nrounded`n'}
		global total_sample_NR = $total_sample_NR + ${nNRrounded`n'}
	}
	
	di $total_sample
	di $total_sample_NR