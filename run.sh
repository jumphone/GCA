NORM_EXP=$1
NET_FILE=$2
OUTPUT=$3 
CPU=$4

SEED=321
PERCENT=0.1
PYTHON="python"
RSCRIPT="Rscript"

$PYTHON step2_BuildIndex.py $NET_FILE $NORM_EXP $OUTPUT\.data $PERCENT
$PYTHON step3_CalZmat.py $OUTPUT\.data $OUTPUT\.data.zmat $CPU
$PYTHON step4_rmSymDist.py $OUTPUT\.data.zmat $OUTPUT\.data.zmat.nosym
$RSCRIPT step5_deMix.R $NORM_EXP $OUTPUT\.data.zmat.nosym $OUTPUT\.data.zmat.nosym.gca_result $CPU $SEED


