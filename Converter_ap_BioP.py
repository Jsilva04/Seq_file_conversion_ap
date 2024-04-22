from Bio import SeqIO
import sys

def get_file_format(filename):
    extension = filename.split('.')[-1]
    if extension.lower() == 'fasta':
        return 'fasta'
    elif extension.lower() == 'nexus':
        return 'nexus'
    elif extension.lower() in ['phylip', 'phy']:
        return 'phylip'
    else:
        return 'invalid'

def convert(input_file, output_file):
    input_format = get_file_format(input_file)
    output_format = get_file_format(output_file)
    SeqIO.convert(input_file, input_format, output_file, output_format, "DNA")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python Converter_ap_BioP.py input_file output_file")
        sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]  
convert(input_file, output_file)
