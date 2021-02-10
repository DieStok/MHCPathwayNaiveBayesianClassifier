Version 1.3.11 of the Short Time-series Expression Miner (STEM)
-------------------------------------------------------------------------------
STEM was developed by Jason Ernst, Dima Patek, and Ziv Bar-Joseph. 
Any questions about the software or bugs found should be emailed to 
Jason Ernst (jernst@cs.cmu.edu).
The latest version of STEM will be available at www.sb.cs.cmu.edu/stem.
Consult the file manual.pdf for comprehensive documentation about STEM.

-------------------------------------------------------------------------------
STEM requires Java 1.4 or later to be installed.  Java can be downloaded from
http://www.java.com.

-------------------------------------------------------------------------------
Once Java is installed, to start STEM in Microsoft Windows double click on stem.cmd.  
Otherwise to start STEM, from a command line while in the stem directory enter the 
command

java -mx1024M  -jar stem.jar 

Append "-d defaults.txt" to the above command to have STEM start with its initial
settings be those specified in the defaults.txt file.  
-------------------------------------------------------------------------------
To start STEM with default settings to analyze the sample data included in the 
stem directory, if in Microsoft Windows double click on the file 
stemGuilleminSample.cmd, otherwise from a command line while in the stem
directory type the command:

java -mx1024M -jar stem.jar -d defaultsGuilleminSample.txt

In the stem directory are three sample data files,
g27_1.txt, g27_2.txt, and vaca.txt.  The data in these files
come from the Stanford Microarray Database entry
http://genome-www5.stanford.edu/cgi-bin/publication/viewPublication.pl?pub_no=193 
which is for the paper
Guillemin, K., Salama, N., Tomplins, L., and Falkow, S.
Cag pathogenicity island-specific responses of gastric
epithelial cells to Helicobacter pylori infection.
Proc Natl Acad Sci. USA, 99 15136-15141.

The data in the g27_1.txt and g27_2.txt files are from the G27 TC1 trial 4 
and G27 TC1 trial 5 time courses respectively and are technical replicates 
of a wildtype experiment.  The data in vaca.txt is a comparison knockout
experiment from the vacA trial 3 time course.
-------------------------------------------------------------------------------
Source code and Javadoc api for the source code are included in the source code 
directory.

-------------------------------------------------------------------------------
As of version 1.3.9 STEM is distributed under a GPL 3.0 license.

STEM makes use of the Piccolo toolkit which is distributed under a BSD license.  
More information about Piccolo can be found at www.cs.umd.edu/hcil/piccolo

STEM makes use of the Batik software which is distributed under an Apache license
More information about Batik can be found at http://xmlgraphics.apache.org/batik/.

STEM makes use of the Gene Ontology and gene annotations provided by Gene Ontology
Consortium members.  More information about the Gene Ontology can be found at
www.geneontology.org.
