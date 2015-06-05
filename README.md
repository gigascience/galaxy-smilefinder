# smilefinder

This repository houses the Galaxy wrappers for the tools described by 
[Guiblet et al., (2015) SmileFinder: a resampling-based approach to evaluate signatures of selection from genome-wide sets of matching allele frequency data in two or more diploid populations. GigaScience 4:1](http://www.gigasciencejournal.com/content/4/1/1).

The [GigaToolShed](http://gigatoolshed.net) is used to host these Galaxy tools.

To clone this repository, [Git Large File Storage](https://git-lfs.github.com/)
has to be installed and made available on your path.

## Folder structure

There is one directory for each tool which has been wrapped for use in Galaxy
and these are contained within the `tools` folder.

## Installation

These Galaxy tools are to be used from within a Galaxy server. They can be
automatically installed via the
[GigaScience Tool Shed](http://gigatoolshed.net). All tools can be manually
installed too. Documentation for automatically and manually installing each 
Galaxy tool can be found within the tools' README file.

## Testing

Functional tests for most of the tools can be found in the tool's XML
configuration file.