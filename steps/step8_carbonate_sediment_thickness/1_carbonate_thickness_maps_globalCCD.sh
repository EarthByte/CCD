#!/bin/zsh

# RDM, January 2022

# Make Atlantic Ocean carbonate thickness maps 

gmt set FONT_ANNOT_PRIMARY 11.5p FONT_LABEL 13p
gmt set MAP_FRAME_TYPE fancy FORMAT_GEO_MAP dddF MAP_TICK_PEN 0.75p MAP_FRAME_PEN 0.75p MAP_TICK_LENGTH_PRIMARY=5p
projection=M
width=18
width_scalebar=`echo "scale=2; $width - 1" | bc `
cpt=carbonate_thickness_blue_orange_red.cpt
# coastlines ------
coastline_col=grey85
sunangle=45

# ffmpeg video settings
rate=6
size=1132x704  # 'yuv420p' requires dimensions to be divisible by 2, so numbers here made even. Otherwise, making the output 50% of the original file
size=2264x1408
quality=25


count=1
for region in global
do
	
	imagedir=carbonate_thickness_maps_$region

	if [ ! -d $imagedir ]; then
		mkdir $imagedir
	fi

# Set region extent based on region
	
frame=-90/90/-180/180


#  loop over reconstructed ages -----------------------------------

age=0
maxage=170
# maxage=1

	while (( ${age} <= $maxage ))
  	do
		echo "Age = $age"

		psfile=$imagedir/compacted_carb_thick_${age}.ps

# input netcdf file with compacted carbonate thickness
		
		grdfile=sediment_thickness_$region/compacted_sediment_thickness_0.25_$age.nc			
		gradfile=sediment_thickness_$region/compacted_sediment_thickness_0.25_${age}_grad.nc

# coastline inputs -----------------------------------
		coastlines=./gplates_export/coastlines_${age}Ma.xy  # antarctica has been removed from this
		
# make gradient grid
		gmt grdgradient $grdfile -G$gradfile -A$sunangle -Ne0.3
		
# plot -------------------------------------------

		gmt psbasemap -R$frame -J${projection}${width} -BWSen+t" " -Bpxa10f5 -Bpya10f5 -Xc -Yc -K -P >$psfile
		gmt grdimage -R -J $grdfile -C$cpt -I$gradfile -Ba0 -K -O >> $psfile
		gmt psxy $coastlines -R -J -G${coastline_col}  -K -O >> $psfile
#		gmt psxy gplates_export/scotia_med_masks/reconstructed_scotia_med_masks_$age.00Ma.xy -R -J -W2 -K -O >> $psfile

# plot scalebar
	# add -Y-7.5 and -DJBC for orthographic SH and or -Y-8.6  and -DJMC for NH 
		gmt psscale -R -J -C$cpt -I0.3 -DJBC+w${width_scalebar}c/0.4c+o0/0.6c+h+ef -Bxa50f10+l"Compacted Carbonate Thickness [m]" -I0.45 -N -O -Y-0.5 -K  >> $psfile	
	
# for SH add -Y-13.5 -X0.7, or NH add -Y13.4 -X0.7

echo ". " | gmt pstext -R0/10/0/10 -JX${width}/9 -F+cTL+f22p,Helvetica,black -D-0.8/0.15 -N -Y2.5 -K -O  >> $psfile # so the images are all the same at the end

echo "$age Ma" | gmt pstext -R -J -F+cTL+f22p,Helvetica,black -D0.5/-10 -N -O >> $psfile

echo "Done making $psfile"

		gmt psconvert $psfile -E300 -TG -P -A0.2c
#		open $psfile

		age=$(($age + 1))

	done
	
# Make movie
	
	input=$imagedir/compacted_carb_thick_%d.png
	reversefile=output/compacted_carb_thick_globalCCD_0-66Ma_$region.mp4
	forwardfile=output/compacted_carb_thick_globalCCD_66-0Ma_$region.mp4

	ffmpeg -r $rate -i $input -q:v 1 -crf $quality -vcodec libx264 -pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2"  -s $size $reversefile
	ffmpeg -i $reversefile -vf reverse $forwardfile
	
	
count=$(($count + 1))
	
done


