#!/bin/bash

#SBATCH --cpus-per-task=1             # 1 core(CPU)
#SBATCH --job-name=getstat   # sensible name for the job
#SBATCH --mem=5G                 # Default memory per CPU is 3GB. Increase memory for this sofware because does not handle multiple cpus
#SBATCH --output=getstat.txt  #output
#SBATCH --constraint=avx2

##note on this: this script will consider any 1010 or so as binary automatically
## french converter from decimal to binary and reverse : https://sebastienguillon.com/test/javascript/convertisseur.html 

if [ -z "$1" ]
then
      echo "Error: No file provided. Terminating."
      exit
else
      FILE=$1
      if [ "$FILE" == "help" ]
      then 
            echo -e "This script is used to extract differents statistics from a vcf file."
            echo -e "How to use:\n"
            echo -e "sbatch get_vcf_stat.slurm [file] [options]\n"
            echo -e "[file]  The path to the file to process\n"
            echo "[options]     A number (in decimal or binary format) indicating the fields to extract."
            echo "Field order is the following: CHROM-POS-END-ID-REF-ALT-QUAL-SVTYPE-GT-FREQ."
            exit
      else
            echo "Processing extraction on file $FILE"
            ##get working directory to pt the results in
            ##DIR=$(echo ${FILE%/*})
      fi
fi

echo "Job number $SLURM_JOB_ID, starting at"
date

echo "Use sacct --format=jobid,jobname,partition,alloccpus,elapsed,state,MaxVMSize,ReqMem,node -j $SLURM_JOB_ID to have more information on the job."


echo "This script is used to extract differents statistics from a vcf file."
echo "Field order is the following: CHROM-POS-END-ID-REF-ALT-QUAL-SVTYPE-GT-FREQ."



#### function to replace the bc command

dec_to_bin() {
    local n=$1
    local bin=""
    while [ "$n" -gt 0 ]; do
        bin=$(( n % 2 ))$bin
        n=$(( n / 2 ))
    done
    echo "${bin:-0}"
}


if [ -z "$2" ]
then
      echo "Error: Unable to understand field wanted because no number provided. Terminating"
      exit
else
    INPUTNUMB=$2
    ##check if the number is specified in binary format or decimal
    
            if echo $INPUTNUMB | grep -qE '^[01]+$'; ##if number is binary
            then
                  BINNUMBER=$INPUTNUMB
            else
                  if (($INPUTNUMB>1023))
                  then 
                        echo "Error: Unable to understand field wanted because number provided is too high. Terminating"
                        exit
                  else
                        ##Tranform the decimal number into binary number
                        ##BINNUMBER=$(echo "obase=2;$INPUTNUMB" | bc) ## where $INPUTNUMB is the number
                        BINNUMBER=$(dec_to_bin "$INPUTNUMB")
                  fi
            fi

fi


##first_digit="${INPUTNUMB:0:1}" ##extract first digit
    ##if [ $first_digit -eq 0 ] || [ $first_digit -eq 1 ]; 

## filling the number with 0

# Get the length of the number
length=${#BINNUMBER}

# Add leading zeros to the number until its length is 10
while [ $length -lt 10 ]; do
  BINNUMBER="0$BINNUMBER"
  length=${#BINNUMBER}
done

##echo $BINNUMBER

##checking the fields wanted

if [ ${BINNUMBER:0:1} -eq 0 ]
then
      CHROMFIELD=""
      HCH=""
elif [ ${BINNUMBER:0:1} -eq 1 ]
then
      CHROMFIELD=" %CHROM"
      HCH="CHROM  "
else
      echo "Error."
      exit
fi

if [ ${BINNUMBER:1:1} -eq 0 ]
then
      POSFIELD=""
      HPOS=""
elif [ ${BINNUMBER:1:1} -eq 1 ]
then
      POSFIELD=" %POS"
      HPOS="POS   "
else
      echo "Error."
      exit
fi

if [ ${BINNUMBER:2:1} -eq 0 ]
then
      ENDFIELD=""
      HEND=""
elif [ ${BINNUMBER:2:1} -eq 1 ]
then
      ENDFIELD=" %END"
      HEND="END   "

else
      echo "Error."
      exit
fi

if [ ${BINNUMBER:3:1} -eq 0 ]
then
      IDFIELD=""
      HID=""
elif [ ${BINNUMBER:3:1} -eq 1 ]
then
      IDFIELD=" %ID"
      HID="ID     "
else
      echo "Error."
      exit
fi

if [ ${BINNUMBER:4:1} -eq 0 ]
then
      REFFIELD=""
      HREF=""
elif [ ${BINNUMBER:4:1} -eq 1 ]
then
      REFFIELD=" %REF"
      HREF="REF   "
else
      echo "Error."
      exit
fi

if [ ${BINNUMBER:5:1} -eq 0 ]
then
      ALTFIELD=""
      HALT=""
elif [ ${BINNUMBER:5:1} -eq 1 ]
then
      ALTFIELD=" %ALT"
      HALT="ALT   "
else
      echo "Error."
      exit
fi

if [ ${BINNUMBER:6:1} -eq 0 ]
then
      QUALFIELD=""
      HQUAL=""
elif [ ${BINNUMBER:6:1} -eq 1 ]
then
      QUALFIELD=" %QUAL"
      HQUAL="QUAL "
else
      echo "Error."
      exit
fi

if [ ${BINNUMBER:7:1} -eq 0 ]
then
      SVTYPEFIELD=""
      HSVTYPE=""
elif [ ${BINNUMBER:7:1} -eq 1 ]
then
      SVTYPEFIELD=" %SVTYPE"
      HSVTYPE="SVTYPE   "
else
      echo "Error."
      exit
fi

if [ ${BINNUMBER:8:1} -eq 0 ]
then
      GTFIELD=""
elif [ ${BINNUMBER:8:1} -eq 1 ]
then
      GTFIELD=" [ %GT]"
else
      echo "Error."
      exit
fi

COMMAND="${CHROMFIELD}${POSFIELD}${ENDFIELD}${IDFIELD}${REFFIELD}${ALTFIELD}${QUALFIELD}${SVTYPEFIELD}${GTFIELD}\n"

echo "Extracting information using the command $COMMAND"

module load BCFtools
bcftools query -f "$COMMAND" $FILE -o ${FILE}_${INPUTNUMB}.stat
## add header if GT not asked
if [ ${BINNUMBER:8:1} -eq 0 ]
then
      HEADER="${HCH}${HPOS}${HEND}${HID}${HREF}${HALT}${HQUAL}${HSVTYPE}"
      sed  -i "1i $HEADER" ${FILE}_${INPUTNUMB}.stat
elif [ ${BINNUMBER:8:1} -eq 1 ]
then
      VAR=""
else
      echo "Error."
      exit
fi

if [ ${BINNUMBER:9:1} -eq 0 ]
then
      VAR=""
elif [ ${BINNUMBER:9:1} -eq 1 ]
then
      echo "Extracting variant frequency"
      module load VCFtools
      vcftools --gzvcf $FILE --freq --out ${FILE}_${INPUTNUMB}.vaf
else
      echo "Error."
      exit
fi

echo "Job $SLURM_JOB_ID ended at"
date
