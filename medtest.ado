/* Function: 有调节的中介效应检验函数         */
/* Date: 23Mar2022                            */
/* Author: Alpha Statistics Studio            */
/* Return: te, ide, me/te                     */

capture program drop medtest
program medtest, rclass
	syntax varlist(max=1), iv(string) mv(string) absorb(string) /* 
						*/ [type(string) cv(string) mo(string) /* 
						*/  mo_value(string) center cluster(string) vce(string)]
	
	if ("`mo'" == ""){
		* direct effect model
		reghdfe `varlist' `iv' `cv', absorb(`absorb') cluster(`cluster') vce(`vce')
		local te = _b[`iv']
		est store m1

		* mediating model
		reghdfe `mv' `iv' `cv', absorb(`absorb') cluster(`cluster') vce(`vce')
		local a = _b[`iv']
		est store m2
		
		* full model
		reghdfe `varlist' `iv' `mv' `cv', absorb(`absorb') cluster(`cluster') vce(`vce')
		local b = _b[`mv']
		local me = `a'*`b'
		local ide = _b[`iv']+`me'
		local degree = `me' / `te'
		est store m3
	}
	if ("`mo'" != ""){
		* generate interaction term
		if ("`center'" != ""){
			qui sum `iv'
			local avg1 = r(mean)
			qui sum `mo'
			local avg2 = r(mean)		
			qui sum `mv'
			local avg3 = r(mean)
			gen `iv'_`mo'=(`iv'-`avg1')*(`mo'-`avg2')
			gen `mv'_`mo'=(`mv'-`avg3')*(`mo'-`avg2')
		}
		if ("`center'" == ""){
			gen `iv'_`mo'=`iv'*`mo'
			gen `mv'_`mo'=`mv'*`mo'
		}
		
		* direct effect model
		reghdfe `varlist' `iv' `mo' `iv'_`mo' `cv', absorb(`absorb') cluster(`cluster') vce(`vce')
		local r1 = _b[`iv']
		local r3 = _b[`iv'_`mo']
		est store m1

		* (X → M) is adjusted by Z
		if ("`type'" == "1"){
			* mediating model
			reghdfe `mv' `iv' `mo' `iv'_`mo' `cv', absorb(`absorb') cluster(`cluster') vce(`vce')
			local r4 = _b[`iv']
			local r6 = _b[`iv'_`mo']
			est store m2
			
			* indirect effect model	with mediating variable
			reghdfe `varlist' `iv' `mv' `mo' `iv'_`mo' `cv', absorb(`absorb') cluster(`cluster') vce(`vce')
			local r7 = _b[`iv']
			local r8 = _b[`mv']
			local r10 = _b[`iv'_`mo']
			est store m3
		}
		
		* (M → Y) is adjusted by Z
		if ("`type'" == "2"){
			* mediating model
			reghdfe `mv' `iv' `cv', absorb(`absorb') cluster(`cluster') vce(`vce')
			local r4 = _b[`iv']
			est store m2
			
			* indirect effect model	with mediating variable
			reghdfe `varlist' `iv' `mv' `mo' `mv'_`mo' `cv', absorb(`absorb') cluster(`cluster') vce(`vce')
			local r7 = _b[`iv']
			local r8 = _b[`mv']
			local r11 = _b[`mv'_`mo']
			est store m3
		}	
		
		* (X → M) & (M → Y) are adjusted by Z	
		if ("`type'" == "3"){
			* mediating model
			reghdfe `mv' `iv' `mo' `iv'_`mo' `cv', absorb(`absorb') cluster(`cluster') vce(`vce')
			local r4 = _b[`iv']
			local r6 = _b[`iv'_`mo']
			est store m2
			
			* indirect effect model	with mediating variable
			reghdfe `varlist' `iv' `mv' `mo' `iv'_`mo' `mv'_`mo' `cv', absorb(`absorb') cluster(`cluster') vce(`vce')
			local r7 = _b[`iv']
			local r8 = _b[`mv']
			local r10 = _b[`iv'_`mo']
			local r11 = _b[`mv'_`mo']
			est store m3
		}


		if ("`mo_value'" != ""){
			* estimation of total effect: {r1+r3*mo}
			local te = `r1'+`r3'*`mo_value'
			
			* estimation of indirect effect:
			if ("`type'" == "1"){
				* {r7+r4*r8+(r6*r8+r10)*mo}
				local ide = `r7'+`r4'*`r8'+(`r6'*`r8'+`r10')*`mo_value'
				local me = (`r4'+`r6'*`mo_value')*(`r8')
			}
			if ("`type'" == "2"){
				* {r7+r4*r8+r4*r11*mo}
				local ide = `r7'+`r4'*`r8'+`r4'*`r11'*`mo_value'
				local me = (`r4')*(`r8'+`r11'*`mo_value')
			}		
			if ("`type'" == "3"){
				* {r7+r4*r8+(r6*r8+r10+r4*r11+r6*r11*mo)*mo}
				local ide = `r7'+`r4'*`r8'+(`r6'*`r8'+`r10'+`r4'*`r11'+`r6'*`r11'*`mo_value')*`mo_value'
				local me = (`r4'+`r6'*`mo_value')*(`r8'+`r11'*`mo_value')			
				local turning = -1*(`r6'*`r8'+`r10'+`r4'*`r11')/(2*`r6'*`r11')
				return scalar turning = round(`turning',.0001)
			}
			
			* mediating degree: ide/te
			local degree = `me' / `te'
			
		}
		drop `iv'_`mo' `mv'_`mo'
	}
	
	di as smcl `"The regression results can be displaied by runing: {stata "esttab m1 m2 m3, compress nogap cells(b(star fmt(%9.4f)) se(par)) stats(N r2_a F,fmt(%9.0g %9.4f))": esttab m1 m2 m3}"'
				
	* return function
	return scalar te = round(`te',.0001)
	return scalar ide = round(`ide',.0001)
	return scalar me = round(`me',.0001)
	return scalar degree = round(`degree',.0001)

end

/*

bootstrap me = r(me) ide = r(ide), rep(100) seed(12345): ///
medtest y11, iv(x12) mv(me2) ///
absorb(pro year) cv(c1 c2 c3 c4 c5 c6 c7) mo(mo1) mo_value(1) type(3) 

esttab  m1 m2 m3, ///
		cells(b(star fmt(%9.4f)) se(par))                                ///
		stats(N r2_a F,fmt(%9.0g %9.4f))     ///
		legend label collabels(none) varlabels(_cons Constant)           ///
		star(* 0.10 ** 0.05 *** 0.01) nogap compress                     ///
		indicate("controls=c1 c2 c3 c4 c5 c6 c7")
*/