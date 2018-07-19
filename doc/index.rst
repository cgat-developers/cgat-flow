.. _contents:

==================================================================
CGATPipelines |version| - NGS and Genomics pipelines
==================================================================

This document brings together the various pipelines and scripts
written before and during the `CGAT Training Programme`_. This
documentation has two parts. Below is the general documentation
covering the complete code collection.

Overview
========

This is the code collecton written by members of the `CGAT Training Programme`_
for the execution of pipelines for next-generation sequencing analysis (NGS).

.. _cgatpipelines:

Documentation
=============

.. toctree::
   :maxdepth: 2

   PipelinesBackground.rst
   InstallingPipelines.rst
   UsingPipelines.rst
   Tutorials.rst
   BuildingPipelines.rst
   PipelineReports.rst
   Reference.rst
   Release.rst
   styleguide.rst
   Developers.rst

For information on how to contribute to the pipeline collection,
please see the `CGAT code collection
<https://www.cgat.org/downloads/public/cgat/documentation/>`_.

Active Pipelines
================

The pipelines in this section are being used and maintained by CGAT. We have travis
`integration <https://travis-ci.org>`_ testing of all of our scripts and we run
full local testing of our pipelines on `jenkins <https://jenkins.io>`_ each evening.

NGS Pipelines
-------------

.. toctree::
   :maxdepth: 1
   :caption: NGS Pipelines
   :name: NGS Pipelines

   readqc: Read QC and processing before mapping <pipelines/pipeline_readqc.rst>
   mapping: Short read mapping <pipelines/pipeline_mapping.rst>
   bamstats: Mapping QC statistics <pipelines/pipeline_bamstats.rst>
   rnaseqdiffexpression: RNA-seq differential expression <pipelines/pipeline_rnaseqdiffexpression.rst>
   peakcalling: ChIP-seq Peak calling <pipelines/pipeline_peakcalling.rst>
   windows: Genomic read distribution <pipelines/pipeline_windows.rst>
   exome: Exome-seq <pipelines/pipeline_exome.rst>
   motifs: motif discovery on intervals <pipelines/pipeline_motifs.rst>

Genomics Pipelines
------------------

.. toctree::
   :maxdepth: 1	

   genesets: Building genomic annotations <pipelines/pipeline_annotations.rst>
   intervals: Annotating genomic intervals <pipelines/pipeline_intervals.rst>

Meta pipelines
--------------

.. toctree::
   :maxdepth: 1	

   testing: Regression testing of pipelines <pipelines/pipeline_testing.rst>

Disclaimer
==========

This collection of pipelines is the outcome of 10 years working in various 
fields in bioinformatics. It contains both the good, the bad and the ugly. 
Use at your own risk.

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

