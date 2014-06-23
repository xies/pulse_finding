#!/bin/bash
set -e

JOB=simulate_pulsing

MAIL_TYPE=e
MAIL_USER=xies@mit.edu

PULSE_DIR='/home/xies/pulse_simulations/pulse_files'
MC_DIR='/home/xies/pulse_simulations/mc_files'

function submit_pulse_simulation {

	INPUT_PULSES=${1}
	ITER_ID=${2}

	OUT_DIR=${MC_DIR}/${INPUT_PULSES}.iter_$ITER_ID
	MAT_OUT=${OUT_DIR}_permutation.mat
	STD_OUT=${OUT_DIR}.stdout
	STD_ERR=${OUT_DIR}.stderr
	SCRIPT=${OUT_DIR}.qsub

	CMD="matlab -nojvm -nodisplay -nodesktop -nosplash -r \"addpath(genpath('~/Code'));batch_simulate_pulsing('${PULSE_DIR}/${INPUT_PULSES}','${MAT_OUT}')\""

	if [ -e ${SCRIPT} ]
	then
		return
	fi

	echo "Writing script to ${SCRIPT}"
	cat << EOF > ${SCRIPT}
#!/bin/bash
#$ -N ${JOB}
#$ -o ${STD_OUT}
#$ -e ${STD_ERR}
#$ -m ${MAIL_TYPE}
#$ -M ${MAIL_USER}

${CMD}
EOF

	echo "Submitting: ${INPUT_PULSES} #$ITER_ID to be simulated..."
	qsub ${SCRIPT}
}

#for (( i=$beg_simulation; i<=$end_simulation; i++ ))
#do
#	submit_pulse_simulation twist.mat $i
#done

submit_pulse_simulation dataA.
