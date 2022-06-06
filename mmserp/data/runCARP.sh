#!/bin/bash 

# =============================================================================
# carp_script.sh
# - Author: Sam Coveney
# - Date: 17 Dec 2020
# -----------------------------------------------------------------------------
# 
# Run mMS simulation with homogeneous parameters as
#   ./carp_script.sh D tau_in tau_out tau_open tau_close S1 PACESITE
#
# Run mMS simulation with heterogeneous parameters as
#   ./carp_script.sh S1 PACESITE
# where the parameters are defined in mesh.par
#
#
# =============================================================================

# pacing is define in file "${PACESITE}.vtx"

#{{{ process bash command line
if [[ $# -eq 2 ]] # hetero simulation
then

    HETERO=1
    echo "Running HETERO-geneous simulation"
    
    # pacing
    S1=${1}
    PACESITE=${2}

elif [[ $# -eq 7 ]] # homo simulation
then

    HETERO=0
    echo "Running HOMO-geneous simulation"

    # pacing
    S1=${6}
    PACESITE=${7}

    # parameters
    DIFFUSION=${1}
    TAU_IN=${2}
    TAU_OUT=${3}
    TAU_OPEN=${4}
    TAU_CLOSE=${5}

else # exit, because command line input is wrong
    echo "[ERROR]: call as <script_name> D t_in t_out t_open t_close S1 PACESITE for homogeneous simulation"
    echo "[ERROR]: call as <script_name> S1 PACESITE for heterogeneous simulation"
    exit
fi
#}}}


#{{{ mesh, param file and output dirs
#====================================

MESH='mesh'
ONAME='./output_'${PACESITE}"_"${S1}
#RSTNAME=${ONAME}'/rst'


# define parameter input
if (( ${HETERO} == 1 ))
then

    # parameters will be defined by parameters file
    PARAMS='+F ./mesh.par'

else

    IMP_GENERAL="-num_imp_regions 1 -imp_region[0].name MYO -imp_region[0].cellSurfVolRatio 1.0 -imp_region[0].im mMS"
    IMP_PARAM="-imp_region[0].im_param V_gate=0.1,a_crit=0.1,tau_in=$TAU_IN,tau_out=$TAU_OUT,tau_open=$TAU_OPEN,tau_close=$TAU_CLOSE"
    G_GENERAL="-num_gregions 1 -gregion[0].name MYO_0 -gregion[0].num_IDs 1 -gregion[0].ID[0] 1"

    # when bidomain equivalent, we need to multiply the diffusion by 2.0
    G_VAL=`python -c "print('{0:0.4f}'.format(2.0*$DIFFUSION*1000.0))"`
    #echo "G_VAL is $G_VAL"

    G_PARAM="-gregion[0].g_el $G_VAL -gregion[0].g_et $G_VAL -gregion[0].g_il $G_VAL -gregion[0].g_it $G_VAL"
    PARAMS="$IMP_GENERAL $IMP_PARAM $G_GENERAL $G_PARAM"

    # specifically name output directory with '_homo
    ONAME=${ONAME}"_homo"

fi

#}}}


#{{{ settings Desktop vs HPC
#===========================

TIMESTEP=100  # time step (microseconds)
ODESUB=5      # nb of sub-iterations for ionic current; can leave 1 now

NPROC=4
CARP='openCARP'
LAUNCHER='mpirun'

CARPCOMMOSTRING="${LAUNCHER} -n ${NPROC} ${CARP} ${PARAMS} -meshname ${MESH} -dt ${TIMESTEP} -ode_fac ${ODESUB}"

#}}}


#{{{ define pacing stimulus
#==========================
# stimulus types: 0=transmembrane; 1=extracell current

# number of beats
NUM_STIM=8

ISTIM=2.0   # intensity
TDUR=1.0    # duration in ms

NPLSS1=1    # nb of stimuli applied to the 1st stimulus
BCLS1=100   # period for stimulus 1

TEND=$((${NUM_STIM}*$S1+300))     # end time (milliseconds) # FIXME: these is S1 + 300 time left after last beat 
SPACEDT=${TEND}     # writes the solution in a file every x milliseconds (5 in this case)

TOTAL_STIMULUS=""

for i in $(seq 0 $((${NUM_STIM}-1))); do

    S1START=$(($S1*$i))

    STIM_LOC="-stimulus[$i].vtx_file ${PACESITE}"

    STIM_COMMON="-stimulus[$i].stimtype 0 -stimulus[$i].strength ${ISTIM} -stimulus[$i].duration ${TDUR} -stimulus[$i].npls ${NPLSS1} -stimulus[$i].bcl ${BCLS1} -stimulus[$i].dump_vtx_file 0 -stimulus[$i].start ${S1START} $STIM_LOC"

    TOTAL_STIMULUS=${TOTAL_STIMULUS}" $STIM_COMMON $STIM_LOC -stimulus[$i].name S${i} " # NOTE: needs whitespace on end

done
#}}}


#{{{ LAT recording
#=================
# all:  0: only 1st activation;
#       1: all the activations (2 columns: node ID and LAT)
# mode: 0 detects depolarisation
#       1 detects repolarisation

TSHACT=0.7  # Threshold to evaluate DEPOLARISATION (i.e. LAT) 
TSHREP=0.1  # Threshold to evaluate REPOLARISATION 90%
TSH50=0.5  # Threshold to evaluate REPOLARISATION 50%
TSH30=0.7  # Threshold to evaluate REPOLARISATION 30%
TSH20=0.8  # Threshold to evaluate REPOLARISATION 20%

LATSSTRING1="-lats[0].measurand 0 -lats[0].all 1 -lats[0].method 1 -lats[0].mode 0 -lats[0].threshold ${TSHACT} -lats[0].ID tact_${TSHACT}"
LATSSTRING2="-lats[1].measurand 0 -lats[1].all 1 -lats[1].method 1 -lats[1].mode 1 -lats[1].threshold ${TSHREP} -lats[1].ID trep_${TSHREP}"
LATSSTRING3="-lats[2].measurand 0 -lats[2].all 1 -lats[2].method 1 -lats[2].mode 1 -lats[2].threshold ${TSH50}  -lats[2].ID trep_${TSH50}"
LATSSTRING4="-lats[3].measurand 0 -lats[3].all 1 -lats[3].method 1 -lats[3].mode 1 -lats[3].threshold ${TSH30}  -lats[3].ID trep_${TSH30}"
LATSSTRING5="-lats[4].measurand 0 -lats[4].all 1 -lats[4].method 1 -lats[4].mode 1 -lats[4].threshold ${TSH20}  -lats[4].ID trep_${TSH20}"
#}}}


#{{{ run CARP simulation
#=======================

# writes screen output every 'TIMEDT' milliseconds
TIMEDT=1    

# three identical S1 beats
${CARPCOMMOSTRING} -simID ${ONAME} -tend ${TEND}.0 -spacedt ${SPACEDT} -timedt ${TIMEDT} -num_stim ${NUM_STIM} ${TOTAL_STIMULUS} -num_LATs 5 ${LATSSTRING1} ${LATSSTRING2} ${LATSSTRING3} ${LATSSTRING4} ${LATSSTRING5}

#}}}


