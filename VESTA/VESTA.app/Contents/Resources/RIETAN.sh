#!/bin/bash
# Execute RIETAN
# The environment variable RIETAN should be defined before execution
COMMANDS=/Applications/RIETAN_VENUS/commands_common
f=${1##*/}
f=${f%.ins}
SAMPLE_DIR=${1%/*}
RIETAN=${2%/*}
export RIETAN
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:.
export LD_LIBRARY_PATH
PATH=$PATH:/usr/local/bin:/usr/texbin
export PATH
cd "$SAMPLE_DIR"
rm $f.itx $f.xyz $f.plt $f.gpd $f.lst $f.pdf
# command open is moved here to show a graph after executing this shell script
open "$SAMPLE_DIR"
$2 $f.ins $f.int $f.bkg $f.itx $f.hkl $f.xyz $f.fos $f.ffe $f.fba $f.ffi $f.ffo $f.vesta $f.plt $f.gpd $f.alb $f.prf $f.inflip $f.exp > $f.lst

if [ -e $f.plt ]; then
   if [ -e /Applications/Gnuplot.app/Contents/Resources/bin/gnuplot-run.sh ]; then
      cat $f.plt | /Applications/Gnuplot.app/Contents/Resources/bin/gnuplot-run.sh > $f.pdf
   else
      gnuplot $f.plt | pstopdf -i -o $f.pdf
      if [ -e /usr/texbin/pdfcrop ]; then
         /usr/texbin/pdfcrop --margins 5 --pdfversion 1.4 $f.pdf $f-crop.pdf
         if [ -e $f-crop.pdf ]; then
            rm $f.pdf
            mv $f-crop.pdf $f.pdf
         fi
      fi
   fi
   if [ -e $f.pdf ]; then
         open -a Preview $f.pdf
   fi
elif [ -e $f.itx ]; then
	if [ -z $3 ] || [ "$3" = "\"" ]
	 then
		open "$SAMPLE_DIR/$f.itx"
		open "$SAMPLE_DIR"
	else
		dir=${3%/*}
		ext=${3##*/}
		if [ "$ext" = "PowderPlot.app" ]
		 then
			#exe=${3##*/}
			cd $dir
			open "$SAMPLE_DIR"
			open $ext
		else
			open -a "$3" "$SAMPLE_DIR/$f.itx"
			open "$SAMPLE_DIR"
		fi
	fi
else
	open -e "$SAMPLE_DIR/$f.lst"
fi
#sleep 5
#rm $SAMPLE_DIR/$f.*
