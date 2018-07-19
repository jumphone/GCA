NORM_EXP=../COMBINE.txt.GBM
NET_FILE=trrust_data/trrust_rawdata.human.tsv
OUTPUT=../COMBINE.txt.GBM.GCA
CPU=4

SEED=123
PERCENT=0.5
PYTHON="python"
RSCRIPT="/home/zhangfeng/bin/bin/Rscript"
CUTOFF=0.3

$PYTHON step1_BuildIndex.py $NET_FILE $NORM_EXP $OUTPUT\.data $PERCENT
$PYTHON step2_CalZmat.py $OUTPUT\.data $OUTPUT\.data.zmat $CPU
$RSCRIPT step3_deMix.R $NORM_EXP $OUTPUT\.data.zmat $OUTPUT\.data.zmat.gca_result $CPU $SEED $CUTOFF
