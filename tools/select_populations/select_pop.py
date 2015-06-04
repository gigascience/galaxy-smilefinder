"""
select_pop.py
For selecting populations for analysis by SmileFinder
Peter Li - GigaScience and BGI-HK
"""

import optparse
import os
import sys
from cStringIO import StringIO


def stop_err(msg):
    sys.stderr.write(msg)
    sys.exit()


def main():
    # Parse command line
    parser = optparse.OptionParser()
    parser.add_option("", "--tool_dir", dest="tool_dir")
    parser.add_option("", "--selection", dest="selection")
    parser.add_option("", "--population1", dest="population1")
    parser.add_option("", "--population2", dest="population2")
    parser.add_option("", "--aggregate", dest="aggregate")

    # Outputs
    parser.add_option("", "--pops", dest="pops")
    opts, args = parser.parse_args()

    # To hold output
    buf = StringIO()

    # Get individual IDs for specified populations
    f = open(opts.tool_dir + '/populations.tab', 'r')
    lines = f.readlines()
    f.close()

    if opts.selection == "specify":
        for i, elem in enumerate(lines):
            tokens = lines[i].split()
            if opts.population1 == tokens[2]:
                buf.write(opts.population1 + "\t" + tokens[0] + "\n")

        for i, elem in enumerate(lines):
            tokens = lines[i].split()
            if opts.population2 == tokens[2]:
                buf.write(opts.population2 + "\t" + tokens[0] + "\n")
    else:
        # List of aggregated individuals are delimited by commas
        agg_individuals = str(opts.aggregate).split(",")
        for i, elem_ind in enumerate(agg_individuals):
            for x, elem_lines in enumerate(lines):
                tokens = lines[x].split()
                if agg_individuals[i] == tokens[2]:
                    buf.write(agg_individuals[i] + "\t" + tokens[0] + "\n")

    # Write results into output
    fh = open(opts.pops, 'w')
    fh.write(buf.getvalue())
    fh.close()

    # Check results in output file
    if os.path.getsize(opts.pops) > 0:
        sys.stdout.write('Status complete')
    else:
        stop_err("Problem with Select Populations process")

if __name__ == "__main__":
    main()

