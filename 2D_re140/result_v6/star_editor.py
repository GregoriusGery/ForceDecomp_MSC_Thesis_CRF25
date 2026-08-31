import os
import glob
import csv

def process_csv_file(input_filepath, output_dir):
    # Read the content of the CSV file
    with open(input_filepath, 'r', newline='') as f:
        reader = csv.reader(f)
        try:
            # Extract headers (first row)
            headers = next(reader)
        except StopIteration:
            print(f"  -> Error: {input_filepath} is empty.")
            return
            
        # Extract data rows
        data_rows = list(reader)

    if not data_rows:
        print(f"  -> Warning: No data rows found in {input_filepath}.")
        return

    # --- NEW LOGIC: Delete 2nd column if total columns > 2 ---
    if len(headers) > 2:
        del headers[1]  # Remove the second column from the header
        
        for row in data_rows:
            if len(row) > 2:
                del row[1]  # Remove the second column from the data rows
    # ---------------------------------------------------------

    # Find the index of the "Physical Time (s)" column
    time_idx = -1
    for i, col in enumerate(headers):
        if col.strip() == "Physical Time (s)":
            time_idx = i
            break
            
    if time_idx == -1:
        print(f"  -> Error: 'Physical Time (s)' column not found in headers of {input_filepath}.")
        return

    # Process and fix "Physical Time (s)" values (multiples of 0.00014)
    multiple = 0.0005
    for i in range(len(data_rows)):
        # Skip empty rows or rows that are missing the time column
        if len(data_rows[i]) <= time_idx:
            continue
            
        try:
            val = float(data_rows[i][time_idx])
        except ValueError:
            # If the value isn't a number (e.g., text), skip it
            continue
            
        # Round to the nearest multiple of 0.00014
        rounded = round(val / multiple) * multiple
        rounded = round(rounded, 5) # Clean up floating-point inaccuracies
        
        # Lookbehind constraint: enforce sequence with the previous row
        if i > 0:
            try:
                prev_val = float(data_rows[i-1][time_idx])
                if rounded <= prev_val:
                    rounded = round(prev_val + multiple, 5)
            except (ValueError, IndexError):
                pass # Skip check if previous row data is missing/invalid
                
        # Lookahead constraint: enforce sequence with the next row
        if i < len(data_rows) - 1:
            try:
                next_val = float(data_rows[i+1][time_idx])
                next_expected = round(round(next_val / multiple) * multiple, 5)
                
                if rounded >= next_expected:
                    rounded = round(next_expected - multiple, 5)
                    # Double-check lookbehind hasn't been broken
                    if i > 0:
                        try:
                            prev_val = float(data_rows[i-1][time_idx])
                            if rounded <= prev_val:
                                rounded = round(prev_val + multiple, 5)
                        except (ValueError, IndexError):
                            pass
            except (ValueError, IndexError):
                pass # Skip check if next row data is missing/invalid
                

        # Update the value in the data row
        data_rows[i][time_idx] = str(rounded)
        
    # Write everything back to the new CSV file in the "edited" directory
    base_filename = os.path.basename(input_filepath)
    output_filepath = os.path.join(output_dir, base_filename)
    
    with open(output_filepath, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(headers)
        writer.writerows(data_rows)

    print(f"  -> Success: Processed and saved.")

if __name__ == "__main__":
    current_directory = os.getcwd()
    output_directory = os.path.join(current_directory, "edited")
    
    # Find all .csv files in the current directory
    csv_files = glob.glob("*.csv")
    
    if not csv_files:
        print("No .csv files found in the current directory.")
    else:
        print(f"Found {len(csv_files)} .csv file(s). Starting batch process...\n")
        
        # Create the "edited" folder if it doesn't already exist
        os.makedirs(output_directory, exist_ok=True)
        
        for file in csv_files:
            print(f"Processing '{file}'...")
            process_csv_file(file, output_directory)
            
        print("\nAll files processed successfully!")