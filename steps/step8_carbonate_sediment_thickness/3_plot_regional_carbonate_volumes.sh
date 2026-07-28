#!/bin/zsh

gmt set MAP_TICK_LENGTH_PRIMARY=5p/2.5p
errorcol=230
count=1


# RDM December 2021
# ------ This script plots regional carbonate volumes through time   ------------------------------------
# ------ comparing the results for regional CCDs versus a global CCD ------------------------------------

# get total Atlantic carbonate volume
# regional CCDs
paste output/carb_vol_north_Natl.txt output/carb_vol_central_Natl.txt output/carb_vol_Satl.txt > output/carb_vol_all.txt
paste output/carb_vol_north_Natl_min.txt output/carb_vol_central_Natl_min.txt output/carb_vol_Satl_min.txt > output/carb_vol_all_min.txt
paste output/carb_vol_north_Natl_max.txt output/carb_vol_central_Natl_max.txt output/carb_vol_Satl_max.txt > output/carb_vol_all_max.txt

gawk '{volume_total = $2 + $6 + $10; print $1, volume_total}' output/carb_vol_all.txt >output/carb_vol_total.txt
gawk '{volume_total = $2 + $6 + $10; print $1, volume_total}' output/carb_vol_all_min.txt >output/carb_vol_total_min.txt
gawk '{volume_total = $2 + $6 + $10; print $1, volume_total}' output/carb_vol_all_max.txt >output/carb_vol_total_max.txt

# global CCD
paste output/carb_vol_globalCCD_north_Natl.txt output/carb_vol_globalCCD_central_Natl.txt output/carb_vol_globalCCD_Satl.txt > output/carb_vol_globalCCD_all.txt
gawk '{volume_total = $2 + $6 + $10; print $1, volume_total}' output/carb_vol_globalCCD_all.txt >output/carb_vol_globalCCD_total.txt

# psfile for carbonate carbon mass (Mt)
psfile=output/atlantic_carbonate_carbon_mass.ps

# psfile for carbonate carbon flux (Mt/yr)
psfile_vcr=output/atlantic_carbonate_carbon_flux.ps

maxvolflux=80
frame=0/66/0/$maxvolflux
scale=X-15/9
gmt psbasemap -R$frame -J$scale -B+t" " -X4 -Y9 -P -K >$psfile_vcr

for region in north_Natl central_Natl Satl total
do

# compute rates of change in carbonate volume (VCR = volume change rate) (skip header row 1) in km^3 (divide by 1e9 as initial results are in m^3)
# Conversion for decompacted carbonate volume change in m^3/my to carbon flux in Mt/yr
# 0.41*0.7*2710*0.12*10^-15
# 1 Mt = 1000000000 kg (1x 10^9)
# 1 myr = 10^6 years
# combine these and thus divide by 1e15
# surface porosity is 59% – this is the same for the entire decompacted section, so solid material is 41% – hence multiply by 0.41
# 70% of this material is carbonate on average – hence multiply by 0.7
# Multiply by sediment grain density of carbonate, which is 2710 kg/m3. 
# To calculate how much C that is we need to consider the molar mass of CaCO3 which is 100.1 g/mol. 
# From that 12 g/mol is C, meaning that there is 12% of C in CaCO3. So to get the weight of C from carbonate 327910 x 0.12.

# regional CCD

# mean
	carb_vol_info1=output/carb_vol_$region.txt
	vcr=output/carb_vcr_$region.txt
	echo "Working on ${carb_vol_info1}, creating $vcr"
	gawk ' NR == 1 {print "# age (Ma), carb vol change (km^3)" }; 
				NR == 2 {age1 = $1; vol1 = $2; next};
			 	{ age2 = $1; vol2=$2; if (age2 != age1) VCR = (vol1-vol2)*0.41*0.7*2710*0.12/1e15; else VCR=0;					
			 	if (VCR != 0) print age1, VCR; age1 = $1; vol1 = $2 }' $carb_vol_info1 >$vcr
# min					
	carb_vol_info1=output/carb_vol_${region}_min.txt
	vcr_min=output/carb_vcr_${region}_min.txt
	echo "Working on ${carb_vol_info1}, creating $vcr_min"
	gawk ' NR == 1 {print "# age (Ma), carb vol change (km^3)" }; 
				NR == 2 {age1 = $1; vol1 = $2; next};
				{ age2 = $1; vol2=$2; if (age2 != age1) VCR = (vol1-vol2)*0.41*0.7*2710*0.12/1e15; else VCR=0;					
				if (VCR != 0) print age1, VCR; age1 = $1; vol1 = $2 }' $carb_vol_info1 >$vcr_min
# max					
	carb_vol_info1=output/carb_vol_${region}_max.txt
	vcr_max=output/carb_vcr_${region}_max.txt
	echo "Working on ${carb_vol_info1}, creating $vcr_max"
	gawk ' NR == 1 {print "# age (Ma), carb vol change (km^3)" }; 
				NR == 2 {age1 = $1; vol1 = $2; next};
				{ age2 = $1; vol2=$2; if (age2 != age1) VCR = (vol1-vol2)*0.41*0.7*2710*0.12/1e15; else VCR=0;					
				if (VCR != 0) print age1, VCR; age1 = $1; vol1 = $2 }' $carb_vol_info1 >$vcr_max
										
# global Boss and Wilkinson CCD - use only mean
					
	carb_vol_info2=output/carb_vol_globalCCD_$region.txt
	vcr_globalCCD=output/carb_vcr_globalCCD_$region.txt
	echo "Working on ${carb_vol_info2}, creating $vcr_globalCCD"
	gawk ' NR == 1 {print "# age (Ma), carb vol change (km^3)" }; 
				NR == 2 {age1 = $1; vol1 = $2; next};
			 	{ age2 = $1; vol2=$2; if (age2 != age1) VCR = (vol1-vol2)*0.41*0.7*2710*0.12/1e15; else VCR=0;	
				if (VCR != 0) print age1, VCR; age1 = $1; vol1 = $2 }' $carb_vol_info2 >$vcr_globalCCD
				
########## Plot regional carbonate volumes in units of km^3 * 10^6
	
# Create error envelope files
	vol_min=output/carb_vol_${region}_min.txt
	vol_max=output/carb_vol_${region}_max.txt
	
	# special case for North Atlantic to remove artefact at 52 Ma
	if [[ $count -eq 1 ]]; then  
		gawk 'NR > 1 {if ($1 != 52) print $1, $2*0.41*0.7*2710*0.12/1e15}' $vol_min  > tmp_min.txt
		gawk 'NR > 1 {if ($1 != 52) print $1, $2*0.41*0.7*2710*0.12/1e15}' $vol_max  > tmp_max.txt	
	else	
		gawk 'NR > 1 {print $1, $2*0.41*0.7*2710*0.12/1e15}' $vol_min  > tmp_min.txt
		gawk 'NR > 1 {print $1, $2*0.41*0.7*2710*0.12/1e15}' $vol_max  > tmp_max.txt
	fi
	
	head -1 tmp_min.txt > aaa
	tail -1 tmp_min.txt > ccc
	sort -n -r tmp_min.txt > ddd
	cat aaa  tmp_max.txt ccc ddd > error_${region}.txt
		
	scale=X-15/9
	frame=0/66/0/4000
	carb_vol_info1=output/carb_vol_$region.txt
	gawk 'NR > 1 {print $1, $2*0.41*0.7*2710*0.12/1e15}' $carb_vol_info1 > tmp1.txt
	
	if [[ $count -eq 1 ]]; then  
		
		gmt psbasemap -R$frame -J$scale -B+t" " -X4 -Y9 -P -K >$psfile
		gmt psxy -R -J error_${region}.txt -t50 -L -G$errorcol -O -K >>$psfile
		gmt psxy -R -J tmp1.txt -W2.5,red -O -K >>$psfile
	
	elif [[ $count -eq 2 ]]; then 
	
		gmt psxy -R$frame  -J error_${region}.txt -t50 -L -G$errorcol -O -K >>$psfile
		gmt psxy -R -J tmp1.txt -W2.5,blue -O -K >>$psfile
		
	elif [[ $count -eq 3 ]]; then 
	
		gmt psxy -R$frame  -J error_${region}.txt -t50 -L -G$errorcol -O -K >>$psfile
		gmt psxy -R -J tmp1.txt -W2.5,cyan -O -K >>$psfile
		
	elif [[ $count -eq 4 ]]; then 
	
		gmt psxy -R$frame -J error_${region}.txt -t50 -L -G$errorcol -O -K >>$psfile
		gmt psxy -R -J tmp1.txt -W3,black -O -K >>$psfile
		
# add global ccd results
		gawk 'NR > 1 {print $1, $2*0.41*0.7*2710*0.12/1e15}' $carb_vol_info2 > tmp2.txt
		gmt psxy -R$frame  -J tmp2.txt -W2.5,magenta  -Bxa10f5+l"Age (Ma)" -By1000f500+l"Carbonate carbon mass (Mt)" -BWeSn -O >>$psfile
		
	else
		echo "Something went wrong"
		
	fi
			
	
########## Plot regional carbonate volume change rate
	
	# Regions:  north_Natl=red, central_Natl=blue, Satl=green, total=black
	
	frame=0/66/0/$maxvolflux
	scale=X-15/9
	gawk 'NR > 1 {print $1, $2}' $vcr > tmp3.txt
	gmt gmtinfo $vcr
	
	# Create error envelope file
	
	gawk 'NR > 1 {print $1, $2}' $vcr_min  > tmp_min.txt
	gawk 'NR > 1 {print $1, $2}' $vcr_max  > tmp_max.txt

	head -1 tmp_min.txt > aaa
	tail -1 tmp_min.txt > ccc
	sort -n -r tmp_min.txt > ddd
	cat aaa  tmp_max.txt ccc ddd > error.txt
	
	echo "Plotting data for $region and count = $count"
		
	if [[ $count -eq 1 ]]; then  
		
		# special case for North Atlantic to remove artefact at 52 Ma
		gawk 'NR > 1 {if ($1 != 52) print $1, $2}' $vcr > tmp3.txt
		gawk 'NR > 1 {if ($1 != 52) print $1, $2}' error.txt > tmp4.txt
		mv tmp4.txt error.txt
		gawk 'NR > 1 {if ($1 != 52) print $1, $2}' error.txt > tmp5.txt
		mv tmp5.txt error.txt
		
		gmt psxy -R$frame -J error.txt -L -G$errorcol -Bxa10f5+l"Age (Ma)" -By10g${maxvolflux}+l"Carbonate carbon flux (Mt/yr)" -BWeSn  -O -K >>$psfile_vcr
		gmt psxy -R -J tmp3.txt -W1.5,red -O -K >>$psfile_vcr
		
	elif [[ $count -eq 2 ]]; then 
		
		gmt psxy -R$frame -J error.txt -L -G$errorcol -O -K >>$psfile_vcr
		gmt psxy -R -J tmp3.txt -W1.5,blue -O -K >>$psfile_vcr
		
	elif [[ $count -eq 3 ]]; then 
		
		gmt psxy -R$frame -J error.txt -L -G$errorcol -O -K >>$psfile_vcr
		gmt psxy -R -J tmp3.txt -W1.5,cyan -O -K >>$psfile_vcr
		
	elif [[ $count -eq 4 ]]; then 
		
		gmt psxy -R$frame -J error.txt -L -G$errorcol -O -K >>$psfile_vcr
		gmt psxy -R -J tmp3.txt -W2.5,black -O -K >>$psfile_vcr
		
		# global CCD
		gawk 'NR > 1 {print $1, $2}' $vcr_globalCCD > tmp4.txt
		gmt psxy -R$frame -J tmp4.txt -W1.5,magenta -O -K >>$psfile_vcr
	else
		echo "Something went wrong"
	fi
	count=$(($count + 1))
	
done

# add total carbonate volume change

gmt psbasemap -R$frame -J$scale -BWeSn -O >>$psfile_vcr

gmt psconvert -Tf -A $psfile_vcr
gmt psconvert -Tf -A $psfile

open output/atlantic_carbonate_carbon_flux.pdf
open output/atlantic_carbonate_carbon_mass.pdf

rm tmp*.txt aaa ccc ddd error*txt


