import hashlib
import os
import csv

def calculate_sha1_checksum(file_path):
    sha1 = hashlib.sha1()
    try:
        with open(file_path, 'rb') as file:
            while chunk := file.read(8192):  
                sha1.update(chunk)
    except Exception as e:
        print(f"Error reading file {file_path}: {e}")
    return sha1.hexdigest()

def process_directory(directory_path, output_csv):
    file_list = []
    for root, _, files in os.walk(directory_path):
        for file in files:
            if file.endswith(('.fasta', '.fa')):
                file_list.append(os.path.join(root, file))
    
    with open(output_csv, mode='w', newline='') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(['File Path', 'SHA-1 Checksum'])  
        
        for file_path in file_list:
            checksum = calculate_sha1_checksum(file_path)
            writer.writerow([file_path, checksum])
            print(f"Processed {file_path}")

if __name__ == '__main__':
    directory_path = '/x/assemblies/'
    output_csv = '/x/checksum_results.csv'
    
    process_directory(directory_path, output_csv)
    print(f"Checksum results have been written to {output_csv}")
