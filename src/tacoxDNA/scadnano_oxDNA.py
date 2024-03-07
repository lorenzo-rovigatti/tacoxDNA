#!/usr/bin/env python3

import scadnano as sc
import argparse
import os
    
def main():
    parser = argparse.ArgumentParser(description="scadnano -> oxDNA converter")
    parser.add_argument("scadnano_file", help="The scadnano input design")
    
    args = parser.parse_args()
    scadnano_file = args.scadnano_file
    
    design = sc.Design.from_scadnano_file(scadnano_file)
    configuration, topology = design.to_oxdna_format()
    
    basename = os.path.basename(scadnano_file)
    topology_file = basename + ".top"
    configuration_file = basename + ".oxdna"

    with open(topology_file, "w") as f:
        f.write(topology)
        
    with open(configuration_file, "w") as f:
        f.write(configuration)


if __name__ == '__main__':
    main()
        