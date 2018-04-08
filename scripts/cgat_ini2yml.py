'''
cgat_ini2yml.py 
=======================================

'''

import sys
import re
import CGATCore.Experiment as E


def main(argv=None):
    """script main.

    parses command line options in sys.argv, unless *argv* is given.
    """

    if argv is None:
        argv = sys.argv

    # setup command line parser
    parser = E.OptionParser(version="%prog version: $Id$",
                            usage=globals()["__doc__"])

    parser.add_option("-n", "--dry-run", dest="dry_run", action="store_true",
                      help="dry run, do not delete any files [%default]")

    parser.set_defaults(dry_run=False)

    # add common options (-h/--help, ...) and parse command line
    (options, args) = E.start(parser, argv=argv)

    indent = 0
    for line in options.stdin:

        if not line.strip():
            options.stdout.write("\n")
            continue
        
        if line.startswith("["):
            section = re.search("\[(.*)\]", line).groups()[0]
            if section == "general":
                indent = 0
                options.stdout.write("\n")
            else:
                indent = 4
                options.stdout.write("{}:\n".format(line.strip()[1:-1]))
                
        elif line.startswith("#"):
            options.stdout.write("{}{}".format(" " * indent, line))
            
        elif "=" in line:
            key, val = re.search("(.*)=(.*)", line).groups()

            if "," in val:
                val = "[{}]".format(val)

            options.stdout.write("{}{}: {}\n".format(" " * indent, key, val))
                
        
    # write footer and output benchmark information.
    E.stop()

if __name__ == "__main__":
    sys.exit(main(sys.argv))
