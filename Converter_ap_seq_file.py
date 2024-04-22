from pathlib import Path
import argparse
import re
import sys
import os


def get_input_file_extension(input_file_dir):
    """
    Obtains the input file's type
    """
    with open(input_file_dir) as input_file:
        first_line = input_file.readline().replace("\n", "")
        if first_line.startswith(">"):
            return "fasta"
        if first_line.lower().startswith("#nexus"):
            return "nexus"
        line_args = first_line.split(" ")
        if len(line_args) == 2:
            if line_args[0].isdigit() and line_args[1].isdigit():
                return "phylip"
        return "invalid"


def get_output_file_extension(output_file_dir):
    """
    Obtains the output file's type
    """
    output_file_array = output_file_dir.split(".")
    extension = output_file_array[len(output_file_array) - 1]
    if extension == "fasta":
        return "fasta"
    elif extension == "nexus":
        return "nexus"
    elif len(output_file_array) == 1 and extension == "phylip":
        return "phylip|phy"
    elif len(output_file_array) > 1 and extension == "phy":
        return "phylip|phy"
    return "invalid"


class Sequence():

    def __init__(self, arguments):

        self.data_dict = {}

        self.input_path = None
        self.output_path = None

        self.input_type = None
        self.output_type = None

        if arguments.input is None and arguments.output is None:
            sys.exit("Must provide an input and output files!")
        if arguments.input is None:
            sys.exit("Must provide an input file!")
        if arguments.output is None:
            sys.exit("Must provide an output file!")

        self.input_path = Path(arguments.input)
        self.output_path = Path(arguments.output)

        if not self.input_path.exists():
            sys.exit("Input file does not exist!")

        path_segments = str(self.output_path).split("\\")
        file_segments = path_segments[len(path_segments) - 1].split(".")
        if not self.output_path.parent.exists() or len(file_segments) == 1:
            self.output_path = None

        self.input_type = get_input_file_extension(self.input_path)
        output_info = get_output_file_extension(arguments.output).split("|")

        self.output_type = output_info[0]
        extension = output_info[len(output_info)-1]

        if self.input_type == "invalid" and self.output_type == "invalid":
            sys.exit("Invalid input and output files!")
        if self.input_type == "invalid":
            sys.exit("Invalid input file!")
        if self.output_type == "invalid":
            sys.exit("Invalid output file!")

        if self.output_path is None:
            self.output_path = Path("output." + extension)

        with open(self.input_path) as open_file:
            sequence_length = None
            if self.input_type == "fasta":
                dna_found = False
                species_name = ""
                dna_sequence = ""
                for line in open_file.readlines():
                    line = line.replace("\n", "")
                    if line.startswith(">"):
                        species_name = line[1:]
                        dna_found = True
                    elif line == "":
                        if sequence_length is None:
                            sequence_length = len(dna_sequence)
                        else:
                            if (sequence_length != len(dna_sequence)):
                                sys.exit("Error: All DNA sequences must be the same length!")
                        self.data_dict.update({species_name: dna_sequence})
                        dna_found = False
                        species_name = ""
                        dna_sequence = ""
                    elif dna_found:
                        dna_sequence += line
            else:
                line_index = 7 if self.input_type == "nexus" else 1
                lines = open_file.readlines()
                while line_index < len(lines):
                    line = lines[line_index].replace("\n", "")
                    if line == ";":
                        break
                    segmented_line = re.split(" ", line)
                    species_name = ""
                    index = 0
                    while index < len(segmented_line) - 1:
                        if segmented_line[index] != "":
                            species_name += segmented_line[index]
                            if segmented_line[index+1] != "":
                                species_name += " "
                        index += 1
                    dna_sequence = segmented_line[len(segmented_line) - 1]
                    if sequence_length is None:
                        sequence_length = len(dna_sequence)
                        else:
                           if (sequence_length != len(dna_sequence)):
                            sys.exit("Error: All DNA sequences must be the same length!")
                    self.data_dict.update({species_name: dna_sequence})
                    line_index += 1

    def get_sequence_lenght(self):
        keys = list(self.data_dict.keys())
        if (len(keys) == 0):
            return 0
        return len(self.data_dict.get(keys[0]))

    def export_data(self):
        keys = self.data_dict.keys()
        match self.output_type:
            case "fasta":
                max_sequence_length = max(len(seq) for seq in self.data_dict.values())
                file_string = ""
                for key in keys:
                    species_name = ">" + key
                    dna_sequence = self.data_dict.get(key)
                    padded_sequence = dna_sequence.ljust(max_sequence_length, '-')
                    file_string += species_name + "\n" + padded_sequence + "\n\n"
                with open(self.output_path, "w") as output_file:
                    output_file.write(file_string)

            case "nexus":
                ntax = len(keys)
                max_sequence_length = max(len(seq) for seq in self.data_dict.values())
                file_string = f"#NEXUS\n\nBEGIN DATA;\nDIMENSIONS NTAX={ntax} NCHAR={max_sequence_length};\nFORMAT DATATYPE=DNA MISSING=N GAP=-;\nMATRIX\n\n"
                for key in keys:
                    species_name = key
                    dna_sequence = self.data_dict.get(key)
                    padded_sequence = dna_sequence.ljust(max_sequence_length, '-')
                    file_string += species_name + "     " + padded_sequence + "\n"
                file_string += ";\n\nEND;\n"
                with open(self.output_path, "w") as output_file:
                    output_file.write(file_string)

            case "phylip":
                ntax = len(keys)
                max_sequence_length = max(len(seq) for seq in self.data_dict.values())
                file_string = f"{ntax} {max_sequence_length}\n"
                for key in keys:
                    species_name = key
                    dna_sequence = self.data_dict.get(key)
                    padded_sequence = dna_sequence.ljust(max_sequence_length, '-')
                    file_string += species_name + "   " + padded_sequence + "\n"
                with open(self.output_path, "w") as output_file:
                    output_file.write(file_string)

        print(f'Data successfully converted and written to: "{self.output_path}"')


parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="Input file")
parser.add_argument("-o", "--output", help="Output file")
args = parser.parse_args()

sequence_holder = Sequence(args)
sequence_holder.export_data()
