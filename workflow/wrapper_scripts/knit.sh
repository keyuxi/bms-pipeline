#!/usr/bin/bash
#SBATCH --job-name="knit"
#SBATCH --time=02:00:00
#SBATCH --output=/home/groups/wjg/kyx/bms_analysis/log/knit-%j.out
#SBATCH --error=/home/groups/wjg/kyx/bms_analysis/log/knit-%j.err
#SBATCH --mail-user=kyx@stanford.edu
#SBATCH --mail-type=FAIL,END,START
#SBATCH --partition=normal,wjg,sfgf,biochem
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=36G

PATH=/home/groups/wjg/resources/software/rstudio-server-2022-07/bin/pandoc:$PATH

# set lib
export R_LIBS_USER="${OAK}/kyx/software/rocker-chromatin/macs2"

# get resources
RESOURCES=`grep -E '^#SBATCH' knit.sh | grep -v 'job-name' | grep -v 'partition' | grep -v 'output'`

# usage:
# sbatch --export=ALL,rmd='01-my_analysis.Rmd' knit.sh

singularity exec /home/groups/wjg/kyx/bms_analysis/Singularity/chromatin-base-0.0.4.sif \
    R --no-save -e "rmarkdown::render('${rmd}', 'html_document')"
