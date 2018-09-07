'''Publication.py - publishing pipeline results
===============================================

Publishing refers to the process of taking all the final results of a
pipeline and putting them somewhere to share with others.

Methods in this module are still in development.

Reference
---------

'''


# import re
# import shutil
# import inspect
# import collections
# import brewer2mpl


# def getPipelineName():
#     '''return the name of the pipeline.

#     The name of the pipeline is deduced by the name of the top-level
#     python script. The pipeline name is the name of the script
#     without any path information and the ``.py`` suffix.

#     Returns
#     -------
#     string

#     '''
#     # use the file attribute of the caller
#     for x in inspect.stack():
#         if x[0].f_globals["__name__"] == "__main__":
#             return os.path.basename(x[0].f_globals['__file__'])[:-3]


# def getPublishDestinations(prefix="", suffix=None):
#     """cgat specific method : return path names of directories
#     for publishing.

#     Arguments
#     ---------
#     prefix : string
#         Prefix to add to output directories.
#     suffix : suffix to add to output directories

#     Returns
#     -------
#     dest_report : string
#          Path for report to export
#     dest_export : string
#          Path for files to export
#     """
#     if not prefix:
#         prefix = PARAMS.get("report_prefix", "default")

#     if prefix == "default":
#         prefix = getPipelineName() + "_"

#     if not suffix:
#         suffix = PARAMS.get("report_suffix", "")

#     dest_report = prefix + "report"
#     dest_export = prefix + "export"

#     if suffix is not None:
#         dest_report += suffix
#         dest_export += suffix

#     return dest_report, dest_export


# def publish_report(prefix="",
#                    patterns=[],
#                    project_id=None,
#                    prefix_project="/ifs/projects",
#                    export_files=None,
#                    suffix=None,
#                    subdirs=False,
#                    ):
#     '''publish report into web directory.

#     Links export directory into web directory.

#     Copies html pages and fudges links to the pages in the
#     export directory.

#     If *prefix* is given, the directories will start with prefix,
#     otherwise, it is looked up from the option ``report_prefix``.
#     If report_prefix is "default", the prefix will be derived
#     from the pipeline name. For example, pipeline_intervals will
#     we copied to ``pipeline_intervals_report``.

#     *patterns* is an optional list of two-element tuples (<pattern>,
#     replacement_string).  Each substitutions will be applied on each
#     file ending in .html.

#     If *project_id* is not given, it will be looked up. This requires
#     that this method is called within a subdirectory of PROJECT_ROOT.

#     *export_files* is a dictionary of files to be exported. The key
#     of the dictionary denotes the targetdirectory within the web
#     directory. The values in the dictionary are the files to be
#     linked to in the direcotry. For example::

#         exportfiles = {
#             "bamfiles" : glob.glob( "*/*.bam" ) + glob.glob( "*/*.bam.bai" ),
#             "bigwigfiles" : glob.glob( "*/*.bw" ),
#             }

#     .. note::
#        This function is cgat specific.

#     '''

#     dest_report, dest_export = getPublishDestinations(prefix, suffix)

#     web_dir = PARAMS["web_dir"]

#     if project_id is None:
#         project_id = getProjectId()

#     src_export = os.path.abspath("export")
#     curdir = os.path.abspath(os.getcwd())

#     # substitute links to export and report
#     base_url = "http://www.cgat.org/downloads/%s" % project_id
#     _patterns = [
#         # redirect export directory
#         (re.compile(src_export),
#          "%(base_url)s/%(dest_export)s" % locals()),
#         # redirect report directory
#         (re.compile(curdir),
#          "%(base_url)s/%(dest_report)s" % locals()),
#         (re.compile('(%s)/_static' %
#                     curdir),
#          "%(base_url)s/%(dest_report)s/_static" % locals())]

#     _patterns.extend(patterns)

#     # add intersphinx mapping - this requires that the name
#     # for the interpshinx redirection (key) corresponds to the
#     # export location with an appended "_report".
#     if CONFIG.has_section("intersphinx"):
#         for key, value in CONFIG.items("intersphinx"):
#             _patterns.append((
#                 re.compile(os.path.abspath(value)),
#                 "%(base_url)s/%(key)s_report" % locals()))
#             # check if the target exists in download location
#             intersphinx_target = os.path.join(
#                 web_dir, key + "_report", "objects.inv")
#             if not os.path.exists(intersphinx_target):
#                 E.warn("intersphinx mapping for '%s' does not exist at %s" %
#                        (key, intersphinx_target))

#     def _link(src, dest):
#         '''create links.

#         Only link to existing targets.
#         '''
#         if os.path.exists(dest):
#             os.remove(dest)

#         if not os.path.exists(src):
#             E.warn("%s does not exist - skipped" % src)
#             return

#         # IMS: check if base path of dest exists. This allows for
#         # prefix to be a nested path structure e.g. project_id/
#         if not os.path.exists(os.path.dirname(os.path.abspath(dest))):
#             E.info('creating directory %s' %
#                    os.path.dirname(os.path.abspath(dest)))
#             os.mkdir(os.path.dirname(os.path.abspath(dest)))

#         os.symlink(os.path.abspath(src), dest)

#     def _copy(src, dest):
#         if os.path.exists(dest):
#             shutil.rmtree(dest)
#         if not os.path.exists(src):
#             E.warn("%s does not exist - skipped" % src)
#             return
#         shutil.copytree(os.path.abspath(src), dest)

#     # publish export dir via symlinking
#     E.info("linking export directory in %s" % dest_export)
#     _link(src_export,
#           os.path.abspath(os.path.join(web_dir, dest_export)))

#     # publish web pages by copying
#     E.info("publishing web pages in %s" %
#            os.path.abspath(os.path.join(web_dir, dest_report)))
#     _copy(os.path.abspath("report/html"),
#           os.path.abspath(os.path.join(web_dir, dest_report)))

#     for root, dirs, files in os.walk(os.path.join(web_dir, dest_report)):
#         for f in files:
#             fn = os.path.join(root, f)
#             if fn.endswith(".html"):
#                 with open(fn) as inf:
#                     data = inf.read()
#                 for rx, repl in _patterns:
#                     data = rx.sub(repl, data)
#                 outf = open(fn, "w")
#                 outf.write(data)
#                 outf.close()

#     if export_files:
#         bigwigs, bams, beds = [], [], []

#         for targetdir, filenames in list(export_files.items()):

#             targetdir = os.path.join(web_dir, targetdir)
#             if not os.path.exists(targetdir):
#                 os.makedirs(targetdir)

#             for src in filenames:
#                 dest = os.path.join(targetdir, os.path.basename(src))
#                 if dest.endswith(".bam"):
#                     bams.append((targetdir, dest))
#                 elif dest.endswith(".bw"):
#                     bigwigs.append((targetdir, dest))
#                 elif dest.endswith(".bed.gz"):
#                     beds.append((targetdir, dest))
#                 dest = os.path.abspath(dest)
#                 if not os.path.exists(dest):
#                     try:
#                         os.symlink(os.path.abspath(src), dest)
#                     except OSError as msg:
#                         E.warn("could not create symlink from %s to %s: %s" %
#                                (os.path.abspath(src), dest, msg))

#         # output ucsc links
#         with open("urls.txt", "w") as outfile:
#             for targetdir, fn in bams:
#                 filename = os.path.basename(fn)
#                 track = filename[:-len(".bam")]
#                 outfile.write(
#                     'track type=bam name="%(track)s" '
#                     'bigDataUrl=http://www.cgat.org/downloads/%(project_id)s/%(targetdir)s/%(filename)s\n'
#                     % locals())

#             for targetdir, fn in bigwigs:
#                 filename = os.path.basename(fn)
#                 track = filename[:-len(".bw")]
#                 outfile.write(
#                     'track type=bigWig name="%(track)s" '
#                     'bigDataUrl=http://www.cgat.org/downloads/%(project_id)s/%(targetdir)s/%(filename)s\n'
#                     % locals())

#             for targetdir, fn in beds:
#                 filename = os.path.basename(fn)
#                 track = filename[:-len(".bed.gz")]
#                 outfile.write(
#                     """http://www.cgat.org/downloads/%(project_id)s/%(targetdir)s/%(filename)s\n""" % locals())

#         E.info("UCSC urls are in urls.txt")

#     E.info(
#         "report has been published at http://www.cgat.org/downloads/%(project_id)s/%(dest_report)s" % locals())


# def publish_tracks(export_files,
#                    prefix="",
#                    project_id=None,
#                    project_name=None,
#                    UCSC_ini=None):
#     '''publish a UCSC Track Hub.

#     This method takes a dictionary of file types associated
#     with files. For each file, a link will be created in
#     the upload directory. The track will be stored under
#     a project name, which will be derived from the location
#     of the working directory.

#     Information about the genome, the upload directory, etc. will be
#     taken from the global configuration dictionary.

#     For example, calling the following code in a pipeline executed
#     in .../proj013/mapping::

#         export_files = {
#             "bamfiles": glob.glob("*/*.bam") + glob.glob("*/*.bam.bai"),
#             "bigwigfiles": glob.glob("*/*.bw"),
#         }
#         publish_tracks(export_files)

#     will create a hub file at
#     :file:`<uploaddir>/OBFUSID/mapping/ucsc.hub`, where
#     OBFUSID is the obfuscated directory entry in the cgat
#     download directory for a particular project.

#     If you want to create group tracks and get them to inherit from a
#     parent, you can supply an filename for a UCSC ini file.  The ini
#     file defines two types of parameters, parents and set_features.
#     Parents define containers with a regex to identify the child
#     tracks. Set_features add additional features to all tracks
#     matching a regex. Parent and set_feature parameters are identified
#     by their respective "parent" or "set_features" prefixes.

#     For example, the following UCSC ini "test.ini" will create a
#     parent multiWig track called "Test" with the UCSC options as
#     defined in the values parameter. The values param must be a comma
#     separated list of key:value pairs which are seperated by a single
#     space. The regex param for parent_test defines the child tracks
#     which will be contained within "Test". The optional colour param
#     defines the colours for the child tracks. Colours are defined
#     using the brewer2mpl python module. Colour parameters must contain
#     the name of the pallete followed by the type of pallette.

#     The ini file below also defines a "set_features" parameter,
#     "bigwigs". Set_feature require a value and regex parameter. In
#     this case, the UCSC options in the values parameter will be added
#     to all tracks matching the ".*bigwig$" regex. As above, the values
#     param must be a comma separated list of key:value pairs which are
#     seperated by a single space. As above, an optional colours
#     parameter can also be given.

#     Note: colour palletes have a maximum number of allowable colours.
#     To see the available palletes and their size, run:
#     >import brewer2mpl
#     >brewer2mpl.print_maps()

#     >cat test.ini
#     #######################
#     #######################

#     [parent_test]
#     FIXME:
#     values=container multiWig,bigDataUrl Test,shortLabel \
#     Test,longLabel Test,type bigWig,viewLimits 0:160,visibility \
#     full,aggregate transparentOverlay,showSubtrackColorOnUi \
#     on,windowingFunction maximum,priority 1.2,configurable \
#     on,autoScale on,dragAndDrop subtracks \

#     regex=.*-Saline-.*bw$

#     colour=Blues,Sequential

#     #######################
#     [set_features_bigwigs]

#     values=configurable on,autoScale on,useScore on,visibility full

#     regex=.*bigwig$

#     colour=Oranges,Sequential
#     #######################
#     #######################

#     Arguments
#     ---------
#     export_files : dict
#         Dictionary mapping filetypes to files.
#     prefix : string
#         will be added to each track.
#     project_id : string
#         The project identifier. If not given, it will be taken from
#         the path of the project directory.
#     project_name : string
#         The project name, typically the project number. If not given,
#         it will be taken from the current directory.

#     '''

#     # the import is located here to avoid cyclical dependencies
#     # between Local.py, Pipeline.py and PipelineUCSC.py
#     import cgatPipelines.PipelineUCSC as PipelineUCSC

#     if not prefix:
#         prefix = PARAMS.get("report_prefix", "")

#     if not UCSC_ini:
#         UCSC_ini = PARAMS.get("ucsc_ini", None)

#     web_dir = PARAMS["web_dir"]
#     if project_id is None:
#         project_id = getProjectId()
#     if project_name is None:
#         project_name = getProjectName()

#     src_export = os.path.abspath("export")
#     dest_report = prefix + "report"
#     dest_export = prefix + "export"

#     hubdir = os.path.join(PARAMS["web_dir"], "ucsc")

#     if not os.path.exists(hubdir):
#         E.info("creating %s" % hubdir)
#         os.mkdir(hubdir)

#     # write the UCSC hub file
#     hubfile = os.path.join(hubdir, "hub.txt")
#     genomesfile = os.path.join(hubdir, "genomes.txt")
#     trackdir = os.path.join(hubdir, PARAMS["genome"])
#     trackfile = os.path.join(hubdir, PARAMS["genome"], "trackDb.txt")
#     trackrelpath = os.path.join(PARAMS["genome"], "trackDb.txt")

#     if os.path.exists(hubfile):
#         with IOTools.openFile(hubfile) as infile:
#             hubdata = PipelineUCSC.readUCSCFile(infile)
#     else:
#         hubdata = [('hub', "cgat-" + project_name),
#                    ('shortLabel', "cgat-" + project_name),
#                    ('longLabel', "Data for cgat project %s" % project_name),
#                    ('genomesFile', "genomes.txt"),
#                    ('email', 'andreas.heger@gmail.com')]

#     E.info("writing to %s" % hubfile)
#     with IOTools.openFile(hubfile, "w") as outfile:
#         PipelineUCSC.writeUCSCFile(outfile, hubdata)

#     # create the genomes.txt file - append to it if necessary.
#     if os.path.exists(genomesfile):
#         with IOTools.openFile(genomesfile) as infile:
#             genomes = PipelineUCSC.readUCSCFile(infile)
#     else:
#         genomes = []

#     if ("genome", PARAMS["genome"]) not in genomes:
#         genomes.append(("genome", PARAMS["genome"]))
#         genomes.append(("trackDb", trackrelpath))

#     E.info("writing to %s" % genomesfile)
#     with IOTools.openFile(genomesfile, "w") as outfile:
#         PipelineUCSC.writeUCSCFile(outfile, genomes)

#     # create the track data
#     if not os.path.exists(trackdir):
#         os.mkdir(trackdir)

#     if os.path.exists(trackfile):
#         E.debug('reading existing tracks from %s' % trackfile)
#         with IOTools.openFile(trackfile) as infile:
#             tracks = PipelineUCSC.readTrackFile(infile)
#     else:
#         tracks = []

#     tracks = collections.OrderedDict(tracks)

#     def getName(name):
#         if name.endswith(".bam"):
#             return "bam", name
#         elif name.endswith(".bw") or name.endswith(".bigwig"):
#             return "bigWig", name
#         elif name.endswith(".bb") or name.endswith(".bigbed"):
#             return "bigBed", name
#         else:
#             return None, None

#     for targetdir, filenames in list(export_files.items()):
#         for src in filenames:
#             dest = os.path.join(trackdir, prefix + os.path.basename(src))
#             dest = os.path.abspath(dest)
#             # create a symlink
#             if not os.path.exists(dest):
#                 try:
#                     os.symlink(os.path.abspath(src), dest)
#                 except OSError as msg:
#                     E.warn("could not create symlink from %s to %s: %s" %
#                            (os.path.abspath(src), dest, msg))
#             ucsctype, trackname = getName(os.path.basename(dest))
#             # ignore invalid types and other files (.bai files, ...)
#             if ucsctype is None:
#                 continue
#             tracks[trackname] = (("bigDataUrl", os.path.basename(dest)),
#                                  ("shortLabel", trackname),
#                                  ("longLabel", trackname),
#                                  ("type", ucsctype))

#     if UCSC_ini:
#         UCSC_PARAMS = loadParameters(UCSC_ini)

#         for param, values in UCSC_PARAMS.items():
#             children = []

#             # find "parent" params
#             if re.match("parent_.*_values", param):
#                 make_group = True
#                 name = param.replace("_values", "")
#                 regex = UCSC_PARAMS[name + "_regex"]

#                 for targetdir, filenames in list(export_files.items()):
#                     for src in filenames:
#                         dest = prefix + os.path.basename(src)
#                         if re.match(regex, dest):
#                             children.append(dest)

#             # find "set features" params
#             elif re.match("set_features_.*_regex", param):
#                 make_group = False
#                 regex = UCSC_PARAMS[param]
#                 name = param.replace("_regex", "")
#                 for targetdir, filenames in list(export_files.items()):
#                     for src in filenames:
#                         dest = prefix + os.path.basename(src)
#                         if re.match(regex, dest):
#                             children.append(dest)
#             else:
#                 continue

#             if name + "_colour" in list(UCSC_PARAMS.keys()):
#                 colour, colour_type = UCSC_PARAMS[name + "_colour"].split(",")
#                 try:
#                     colours = brewer2mpl.get_map(
#                         colour, colour_type, max(3, len(children)))
#                 except ValueError as error:
#                     print(("Could not set colours for %s. See error message"
#                            "%s" % (",".join(children), error)))
#                 colours = colours.colors
#                 # make the colours a shade darker
#                 colours = [[max(0, y - 25) for y in x] for x in colours]
#             else:
#                 colours = None

#             for n, child in enumerate(children):
#                 if make_group:
#                     # make a parent and a copy of the child so we have
#                     # two tracks, one grouped, one by itself
#                     values = UCSC_PARAMS[param]
#                     tracks[name] = [x.split(" ") for x in values.split(",")]
#                     group_trackname = child + "_grouped"
#                     tracks[group_trackname] = tracks[child]
#                     tracks[group_trackname] += (("parent", name),)

#                 else:
#                     # just add the values to the child
#                     values = UCSC_PARAMS[name + "_values"]
#                     tracks[child] += tuple([x.split(" ") for x in values.split(",")])

#                 if colours:
#                     rgb = ",".join(map(str, colours[n]))
#                     tracks[child] += (("color", rgb),)
#                     if make_group:
#                         tracks[group_trackname] += (("color", rgb),)

#     E.info("writing to %s" % trackfile)
#     with IOTools.openFile(trackfile, "w") as outfile:
#         PipelineUCSC.writeTrackFile(outfile, list(tracks.items()))

#     E.info(
#         "data hub has been created at http://www.cgat.org/downloads/%(project_id)s/ucsc/hub.txt" % locals())

# def publish_notebooks():
#     '''publish report into web directory.'''

#     dirs = getProjectDirectories()

#     notebookdir = dirs['notebookdir']
#     exportdir = dirs['exportdir']
#     exportnotebookdir = os.path.join(exportdir, "notebooks")

#     if not os.path.exists(exportnotebookdir):
#         os.makedirs(exportnotebookdir)

#     statement = '''
#     cd %(exportnotebookdir)s;
#     ipython nbconvert
#     %(notebookdir)s/*.ipynb
#     --to html
#     ''' % locals()

#     E.run(statement)
