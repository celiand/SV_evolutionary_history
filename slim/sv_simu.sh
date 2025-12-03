#!/bin/bash

#SBATCH --cpus-per-task=2            # 1 core(CPU)
#SBATCH --job-name=svslimsim  # sensible name for the job
#SBATCH --mem=20G                 # Default memory per CPU is 3GB. Increase memory for this sofware because does not handle multiple cpus
#SBATCH --output=svslimsim.txt  #output


WD=Data/slimfolder/


echo "Use sacct --format=jobid,jobname,partition,alloccpus,elapsed,state,MaxVMSize,ReqMem,node -j $SLURM_JOB_ID to have more information on the job."


eval "$(conda shell.bash hook)"


conda activate $WD/slimconda 

cd build

failed_message="simulation failed"
success_message="simulation succeed"

casestudied="recurrent"
selectiontype="selec"



coeffselec=1.1
recurrency_gen=1000 ### one SV occurs every X generation on average.


echo -e "nbrun\tType_of_Test\tSelection\tTypeOfSelection\tIssuccess\tRec_rate">$WD/resultsimu.txt


for ((i=0; i<10000; i++))
    do
        if [ "$casestudied" == "ancient_old" ]; then
            if [ "$selectiontype" == "no_selec" ]; then
                OUTPUTSLIM=$(./slim -d nbrun=$i ../scripts/simulation_ancient_old_sv_no_selec.slim)
                coeffselecprint=0.0
            elif [ "$selectiontype" == "selec_het" ]; then
                OUTPUTSLIM=$(./slim -d nbrun=$i -d coefsel=$coeffselec ../scripts/simulation_ancient_old_sv_selec_het.slim)
                coeffselecprint=$coeffselec
            else
                OUTPUTSLIM=$(./slim -d nbrun=$i -d coefsel=$coeffselec ../scripts/simulation_ancient_old_sv_selec.slim)
                coeffselecprint=$coeffselec
            fi
        elif [ "$casestudied" == "recurrent" ]; then
            if [ "$selectiontype" == "no_selec" ]; then
                OUTPUTSLIM=$(./slim -d nbrun=$i -d rec_rate=$recurrency_gen ../scripts/simulation_reccurent_sv_no_selec.slim)
                coeffselecprint=0.0
            elif [ "$selectiontype" == "selec_het" ]; then
                OUTPUTSLIM=$(./slim -d nbrun=$i -d rec_rate=$recurrency_gen -d coefsel=$coeffselec ../scripts/simulation_reccurent_sv_selec_het.slim)
                coeffselecprint=$coeffselec
            else
                OUTPUTSLIM=$(./slim -d nbrun=$i -d coefsel=$coeffselec -d rec_rate=$recurrency_gen ../scripts/simulation_reccurent_sv_selec.slim)
                coeffselecprint=$coeffselec
            fi
        fi

        message_part="${OUTPUTSLIM#*MESSAGE: }"
    
        ## cut the extra string from the message
        message_part="${message_part%\"}"

        ### check if the SV is still there (success) or if it disappeared during the simulation (failed)
        if [[ "$message_part" == "$failed_message" ]]; then
            ### failed attempts
            echo -e "$i\t$casestudied\t$coeffselecprint\t$selectiontype\tfailed\t$recurrency_gen">>$WD/resultsimu.txt
        elif [[ "$message_part" == "$success_message" ]]; then
            echo -e "$i\t$casestudied\t$coeffselecprint\t$selectiontype\tsuccess\t$recurrency_gen">>$WD/resultsimu.txt
        else
            var="nothing"
        fi

        

    done