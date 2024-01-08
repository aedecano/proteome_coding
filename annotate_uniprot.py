import argparse
import requests
from Bio import SeqIO

def map_fasta_to_uniprot(input_file, output_file):
    print(f"Mapping protein sequences from {input_file} to UniProt IDs...")

    # Read the protein sequences from the FASTA file
    sequences = SeqIO.to_dict(SeqIO.parse(input_file, "fasta"))

    # Initialize a dictionary to store UniProt mappings
    uniprot_mappings = {}

    # UniProt API endpoint for mapping
    uniprot_api_url = "https://www.ebi.ac.uk/proteins/api/coordinates/"

    for seq_id, sequence in sequences.items():
        # Construct the UniProt API request URL for the given sequence
        api_request_url = f"{uniprot_api_url}{seq_id}"

        # Send a GET request to the UniProt API
        response = requests.get(api_request_url)

        # Check if the request was successful (status code 200)
        if response.status_code == 200:
            # Parse the JSON response
            result = response.json()

            # Extract UniProt ID from the response
            uniprot_id = result.get("id", None)

            if uniprot_id:
                # Store the mapping in the dictionary
                uniprot_mappings[seq_id] = uniprot_id
            else:
                print(f"No UniProt ID found for sequence {seq_id}")

    # Write the mapping results to the output file
    with open(output_file, "w") as output_file:
        for seq_id, uniprot_id in uniprot_mappings.items():
            output_file.write(f"{seq_id}\t{uniprot_id}\n")

    print(f"Mapping completed. Results written to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Map protein sequences from a FASTA file to UniProt IDs.")
    parser.add_argument("input_file", help="Path to the input FASTA file")
    parser.add_argument("output_file", help="Path to the output file for UniProt mappings")

    args = parser.parse_args()

    # Call the function with command-line arguments
    map_fasta_to_uniprot(args.input_file, args.output_file)
