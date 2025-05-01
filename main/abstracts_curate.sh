#!/bin/bash
SCRIPT="abstracts_curate.py"
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
LOG_DIR=$(echo "$SCRIPT_DIR" | sed 's/main/log/')

echo "This script runs $SCRIPT in the background."
echo "Make sure you have the right version of JSON Schema and TXT Query in the data_input folder."
echo "----------------------------------------------------------------------------"

DATE=$(date '+%Y%m%d')
LOG="abstracts_curate_${DATE}.log"
LOG_PATH="$LOG_DIR/$LOG"

read -e -p "Enter input directory name:" -i " pubmed_abstracts" INPUT
read -e -p "Enter output directory name:" -i " extracted_abstracts_$DATE" OUTPUT

echo "----------------------------------------------------------------------------"

echo "Job started at: $(date '+%Y-%m-%d %H:%M:%S')."
echo "Running ${SCRIPT} with input: ${INPUT} and output: ${OUPUT}."
nohup python3 -u ${SCRIPT} -i ${INPUT} -o ${OUTPUT} > $LOG_PATH 2>&1 &
PID=$(ps ux | grep ${SCRIPT} | grep python3 | awk '{print $2}')
echo -e "Background process initiated. Your PID is: \e[32m${PID}\e[0m."
echo "Check $LOG to monitor progress. To terminate job, do kill [PID]"
echo -e "To track API spendings, go to OpenAI's [\e[34m\e]8;;https://platform.openai.com/usage\aUsage Dashboard\e]8;;\a\e[0m]"