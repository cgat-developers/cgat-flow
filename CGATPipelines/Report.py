'''Report.py - building reports
===============================

Reference
---------

'''
import os
import sys
import ruffus

import CGATCore.Experiment as E
import CGATCore.Pipeline as P
import CGATCore.IOTools as IOTools


def run_report(clean=True,
               with_pipeline_status=True,
               pipeline_status_format="svg"):
    '''run CGATreport.

    This will also run ruffus to create an svg image of the pipeline
    status unless *with_pipeline_status* is set to False. The image
    will be saved into the export directory.

    '''

    params = P.get_params()

    if with_pipeline_status:
        targetdir = params["exportdir"]
        if not os.path.exists(targetdir):
            os.mkdir(targetdir)

        ruffus.pipeline_printout_graph(
            os.path.join(
                targetdir,
                "pipeline.%s" % pipeline_status_format),
            pipeline_status_format,
            ["full"],
            checksum_level=params["ruffus_checksums_level"]
        )

    dirname, basename = os.path.split(P.get_caller().__file__)

    report_engine = params.get("report_engine", "cgatreport")
    assert report_engine in ('sphinxreport', 'cgatreport')

    docdir = os.path.join(dirname, "pipeline_docs", IOTools.snip(basename, ".py"))
    themedir = os.path.join(dirname, "pipeline_docs", "themes")
    relpath = os.path.relpath(docdir)
    trackerdir = os.path.join(docdir, "trackers")

    # use a fake X display in order to avoid windows popping up
    # from R plots.
    xvfb_command = IOTools.which("xvfb-run")

    # permit multiple servers using -a option
    if xvfb_command:
        xvfb_command += " -a "
    else:
        xvfb_command = ""

    # if there is no DISPLAY variable set, xvfb runs, but
    # exits with error when killing process. Thus, ignore return
    # value.
    # print os.getenv("DISPLAY"), "command=", xvfb_command
    if not os.getenv("DISPLAY"):
        erase_return = "|| true"
    else:
        erase_return = ""

    if os.path.exists("conf.py"):
        conf_dir = os.path.abspath(".")
    else:
        conf_dir = os.path.join(os.path.dirname(__file__), "configuration")

    # in the current version, xvfb always returns with an error, thus
    # ignore these.
    erase_return = "|| true"

    if clean:
        clean = "rm -rf report _cache _static;"
    else:
        clean = ""

    # with sphinx >1.3.1 the PYTHONPATH needs to be set explicitely as
    # the virtual environment seems to be stripped. It is thus set to
    # the contents of the current sys.path
    syspath = ":".join(sys.path)

    statement = '''
    %(clean)s
    (export SPHINX_DOCSDIR=%(docdir)s;
    export SPHINX_THEMEDIR=%(themedir)s;
    export PYTHONPATH=%(syspath)s;
    %(xvfb_command)s
    %(report_engine)s-build
    --num-jobs=%(report_threads)s
    sphinx-build
    -b html
    -d %(report_doctrees)s
    -c %(conf_dir)s
    -j %(report_threads)s
    %(docdir)s %(report_html)s
    >& report.log %(erase_return)s )
    '''

    P.run(statement)

    E.info('the report is available at %s' % os.path.abspath(
        os.path.join(params['report_html'], "contents.html")))

