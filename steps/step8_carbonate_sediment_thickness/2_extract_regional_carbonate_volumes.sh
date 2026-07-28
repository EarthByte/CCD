#!/bin/zsh

# RDM December 2021
# ------ This script extracts carbonate volumes and mean thicknesses from carbonate thickness grids 
# ------ and writes them to a simple text file -----------------------------------------------------
# Grid masking is also performed for areas where manual masking is required (regions 2 and 5)
	
mkdir -p output
count=1

# for region in north_Natl
for region in north_Natl central_Natl Satl
do

############### alternatively comment/uncomment two lines below for Boss and Wilkinson CCD or regional CCD ####################
# regional CCD
	dir=sediment_thickness_$region 
	carb_vol_info=output/carb_vol_$region.txt
	carb_dep_area_info=output/carb_dep_area_$region.txt

# Boss and Wilkinson CCD 
#	dir=sediment_thickness_globalCCD_$region 
# 	carb_vol_info=output/carb_vol_globalCCD_$region.txt

# remove the carbinfo file in case it exists, then create an empty file

	if [[ -s $carb_vol_info ]]; then
		rm $carb_vol_info
	fi
	
	touch $carb_vol_info
	
	if [[ -s $carb_dep_area_info ]]; then
		rm $carb_dep_area_info
	fi
	
	touch $carb_dep_area_info
	
	
	echo "# Age (Ma), decompacted carbonate sediment volume (m^3), mean decompacted carbonate sediment thickness (m), area (km^2)" >$carb_vol_info
	
	echo "# Age (Ma), carbonate deposition area (km2)" >$carb_dep_area_info
	
	age=0
	maxage=67

	while (( ${age} <= $maxage ))
 	 	do

# ---- Input and output files -----

		grid1=$dir/decompacted_sediment_thickness_0.2_$age.nc
		grid2=$dir/compacted_sediment_thickness_0.2_$age.nc
		
		deposition_mask_grid_raw=$dir/deposition_mask_0.2_$age.nc
		deposition_mask_grid_clean=$dir/deposition_mask_0.2_clean_$age.nc
		echo "Working on $grid1"
		
# 
# apply masks to remove artefacts: Gulf of Mexico and Caribbean, Mediterranean, Scotia Sea
# apply to both decompacted and compacted carbonate grids
		
		if [[ $count -eq 1 ]]; then  
			
			frame=-78/18/44/78
			maskfile=gplates_export/scotia_med_masks/reconstructed_scotia_med_masks_$age.00Ma.xy
			maskgrid=$dir/mask_scotia_med_$age.nc
			gmt grdmask $maskfile -R$frame -I12m -N1/NaN/NaN -G$maskgrid
			grid_masked=$dir/decompacted_sediment_thickness_0.2_masked_$age.nc
			gmt grdmath $grid1 $maskgrid MUL = $grid_masked
			grid_masked=$dir/compacted_sediment_thickness_0.2_masked_$age.nc
			gmt grdmath $grid2 $maskgrid MUL = $grid_masked
			
			# Fix deposition mask grid
			gmt grdmath $deposition_mask_grid_raw $maskgrid MUL = $deposition_mask_grid_clean
	
		elif [[ $count -eq 2 ]]; then  
			
			frame=-82/14/-10/44
			maskfile1=gplates_export/gom_carib_masks/reconstructed_gom_carib_masks_$age.00Ma.xy
			maskfile2=gplates_export/scotia_med_masks/reconstructed_scotia_med_masks_$age.00Ma.xy
			maskgrid=$dir/mask_gom_carib_med_$age.nc
			gmt grdmask $maskfile1 $maskfile2 -R$frame -I12m -N1/NaN/NaN -G$maskgrid
			grid_masked=$dir/decompacted_sediment_thickness_0.2_masked_$age.nc
			gmt grdmath $grid1 $maskgrid MUL = $grid_masked
			grid_masked=$dir/compacted_sediment_thickness_0.2_masked_$age.nc
			gmt grdmath $grid2 $maskgrid MUL = $grid_masked
			
			# Fix deposition mask grid
			gmt grdmath $deposition_mask_grid_raw $maskgrid MUL = $deposition_mask_grid_clean
			
		elif [[ $count -eq 3 ]]; then
			
			frame=-58/20/-56/-10
			maskfile=gplates_export/scotia_med_masks/reconstructed_scotia_med_masks_$age.00Ma.xy
			maskgrid=$dir/mask_scotia_$age.nc
			gmt grdmask $maskfile -R$frame -I12m -N1/NaN/NaN -G$maskgrid
			grid_masked=$dir/decompacted_sediment_thickness_0.2_masked_$age.nc
			gmt grdmath $grid1 $maskgrid MUL = $grid_masked
			grid_masked=$dir/compacted_sediment_thickness_0.2_masked_$age.nc
			gmt grdmath $grid2 $maskgrid MUL = $grid_masked
			
			# Fix deposition mask grid
			gmt grdmath $deposition_mask_grid_raw $maskgrid MUL = $deposition_mask_grid_clean
		fi
		
# write vol in m^3 and mean thickness in meters, use decompacted grids

			gmt grdvolume -Se $dir/decompacted_sediment_thickness_0.2_masked_$age.nc  | awk '{print age, $3, $4, $2/1000000}' age=$age >>$carb_vol_info
			
# write carbonate deposition area

			gmt grdvolume -Sk $deposition_mask_grid_clean  | awk '{printf "%2d\t %10d\n", $1, $2/1e6}' age=$age >>$carb_dep_area_info
			
		age=$(($age + 1))
	done
	
	count=$(($count + 1))
done
