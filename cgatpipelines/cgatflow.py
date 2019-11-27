'''
cgatflow.py - Computational Genomics Analysis Workflows
=======================================================

:Tags: Genomics

To use a specific workflow, type::

    cgatflow <workflow> [workflow options] [workflow arguments]

For this message and a list of available keywords type::

    cgatflow --help

To get help for a specify workflow, type::

    cgatflow <workflow> --help
'''

import os
import sys
import re
import glob
import imp
import cgatpipelines
import collections
import shlex
import subprocess
import cgatcore.iotools as IOTools


def mapKeyword2Script(path):
    '''collect keywords from scripts.'''

    map_keyword2script = collections.defaultdict(list)

    for script in glob.glob(os.path.join(path, "*.py")):
        s = os.path.basename(script)[:-3]
        with IOTools.open_file(script, 'r') as inf:
            data = [x for x in inf.readlines(10000) if x.startswith(':Tags:')]
            if data:
                keywords = [x.strip() for x in data[0][6:].split(' ')]
                for x in keywords:
                    if x:
                        map_keyword2script[x].append(s)

    return map_keyword2script

def printListInColumns(l, ncolumns):
    '''output list *l* in *ncolumns*.'''
    ll = len(l)

    if ll == 0:
        return

    max_width = max([len(x) for x in l]) + 3
    n = ll // ncolumns
    if ll % 3 != 0:
        n += 1

    # build columns
    columns = [l[x * n:x * n + n] for x in range(ncolumns)]

    # add empty fields for missing columns in last row
    for x in range(ncolumns - (len(l) % ncolumns)):
        columns[-(x + 1)].append('')

    # convert to rows
    rows = list(zip(*columns))

    # build pattern for a row
    p = '%-' + str(max_width) + 's'
    pattern = ' '.join([p for x in range(ncolumns)])

    # put it all together
    return '\n'.join([pattern % row for row in rows])


def main(argv=None):

    argv = sys.argv

    if argv[1] == "R":
        path = os.path.join(os.path.abspath(os.path.dirname(cgatpipelines.__file__)),
                            "Rtools")

        if len(argv) == 2 or argv[2] == "--help" or argv[2] == "-h":
            print((globals()["__doc__"]))

            map_keyword2script = mapKeyword2Script(path)

            if len(argv) <= 2:

                print('CGAT R-tools are grouped by keywords. The following keywords')
                print('are defined:\n')
                print(("%s\n" % printListInColumns(list(map_keyword2script.keys()),
                                                   3)))

            elif 'all' in argv[2:]:
                print("The list of all available commands is:\n")
                print(("%s\n" % printListInColumns(
                    sorted([os.path.basename(x)[:-3]
                            for x in glob.glob(os.path.join(path, "*.py"))]),
                    3)))

            else:
                for arg in argv[2:]:
                    if arg in map_keyword2script:
                        print(("Tools matching the keyword '%s':\n" % arg))
                        print(('%s\n' % printListInColumns(
                            sorted(map_keyword2script[arg]),
                            3)))
            return

        command = argv[2]

        command = re.sub("-", "_", command)
        if os.path.exists(os.path.join(path, command + ".R")):
            rscriptname = os.path.join(path, command + ".R")
            statement = (
                "export R_ROOT={r_root} && "
                "Rscript {rscriptname} "
                "{args}".format(
                    r_root=os.path.dirname(path),
                    rscriptname=rscriptname,
                    args=" ".join([shlex.quote(x) for x in argv[3:]])))
            return subprocess.call(statement, shell=True, executable=os.environ["SHELL"])
        else:
        # remove 'cgatflow' and "R" from sys.argv
            del sys.argv[0:2]

            (file, pathname, description) = imp.find_module(command, [path, ])
            module = imp.load_module(command, file, pathname, description)
        module.main(sys.argv)

    path = os.path.join(os.path.abspath(os.path.dirname(cgatpipelines.__file__)), "tools")
    relpath = os.path.abspath("../src")

    paths = [path, relpath]

    if len(argv) == 1 or argv[1] == "--help" or argv[1] == "-h":
        pipelines = []
        for path in paths:
            pipelines.extend(glob.glob(os.path.join(path, "pipeline_*.py")))
        print((globals()["__doc__"]))
        print("The list of all available pipelines is:\n")
        print("{}\n".format(
            printListInColumns(
                sorted([os.path.basename(x)[len("pipeline_"):-len(".py")] for x in pipelines]),
                3)))
        return

    command = argv[1]
    command = re.sub("-", "_", command)
    pipeline = "pipeline_{}".format(command)

    # remove 'cgatflow' from sys.argv
    del sys.argv[0]

    (file, pathname, description) = imp.find_module(pipeline, paths)

    module = imp.load_module(pipeline, file, pathname, description)

    module.main(sys.argv)


if __name__ == "__main__":
    sys.exit(main())
