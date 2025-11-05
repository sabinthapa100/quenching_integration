#!/bin/sh

# List of qhat0 values to run
QHAT_VALUES="0.051 0.09"

echo "Running Quenching for Multiple qhat0 Values"

for val in $QHAT_VALUES; do
  echo "--- Starting run for qhat0 = $val ---"

  # Define the final directory name and a temporary log file
  DIR_NAME="output_qhat0_${val}"
  TEMP_LOG="temp_log_${val}.log"

  # Run the command and use 'tee' to show output AND save to a temp file
  {
    ./quenching -qhat0 "$val"
    echo "Run for qhat0=$val Done!"
  } 2>&1 | tee "$TEMP_LOG"

  # Rename the 'output' directory created by the program
  mv output "$DIR_NAME"

  # Move the temporary log into the newly renamed directory
  mv "$TEMP_LOG" "$DIR_NAME/quenching.log"

  echo "Log file saved to $DIR_NAME/quenching.log"
  echo "--- Finished run for qhat0 = $val ---"
  echo "" # Adds a blank line for better readability
done

echo "All quenching runs are complete."
