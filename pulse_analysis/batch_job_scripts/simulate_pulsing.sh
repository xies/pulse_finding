#!/bin/bash
set -e

JOB=simulate_pulsing

MEMORY=16000
MAIL_TYPE=END
MAIL_USER=xies@mit.edu

DATA_DIR = '/home/xies/pulse_simulation/wt_pulses/'

function submit_pulse_simulation {
	ORIGINAL_PULSES=${1}

	OUT_DIR=${DATA_DIR}/${ORIGINAL_PULSE/}
	STD_OUT=${DATA_DIR}/${ORIGINAL_PULSE}.stdout
	STD_ERR=${DATA_DIR}/${ORIGINAL_PULSE}.stderr
	SCRIPT=${DATA_DIR}/${ORIGINAL_PULSE}.qsub

	CMD = 'matlab -nojvm -nodisplay -nodesktop -nosplash -r "batch_simulate_pulsing(${ORIGINAL_PULSES}); exit;"'

	if [ -e ${SCRIPT} ]
	then
		return
	fi

	cat << EOF > ${SCRIPT}
#!/bin/bash
#$ -N = ${JOB}
#$ -o = ${STD_OUT}
#$ -e = ${STD_ERR}
#$ -m = ${MAIL_TYPE}
#$ -M = ${MAIL_USER}

${CMD}
EOF

	echo "Submitting: ${ORIGINAL_PULSE}"
	sbatch ${SCRIPT}
}

submit_pulse_simulation
