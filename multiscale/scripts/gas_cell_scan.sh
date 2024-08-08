#!/bin/bash

# Directory containing the .h5 files
SOURCE_DIR="."

# Script to run on each file
PROCESS_SCRIPT="$MULTISCALE_SCRIPTS/run_multiscale_SUNRISE.sh"

# Check if the process script exists
if [ ! -f "$PROCESS_SCRIPT" ]; then
    echo "Process script not found: $PROCESS_SCRIPT"
    exit 1
fi

# Loop through each .h5 file in the source directory
for FILE in ${SOURCE_DIR}/results_pressure_*.h5; 
do
    # Get the base name of the file (e.g., results_pressure_1.h5)
    BASE_NAME=$(basename "$FILE")
    
    # Extract the simulation number from the filename using sed
    SIMULATION_NUMBER=$(echo "$BASE_NAME" | sed -n 's/^results_pressure_\([0-9]*\)\.h5$/\1/p')
    
    # Create a directory with the extracted simulation number
    DIR_NAME="t${SIMULATION_NUMBER}"
    mkdir -p "$DIR_NAME"
    
    # Copy the file to the new directory
    cp "$FILE" "$DIR_NAME/"
    
    # Run the process script on the copied file
    (cd "$DIR_NAME" && bash "$PROCESS_SCRIPT")
done

echo "All files processed."
