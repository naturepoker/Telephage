#!/usr/bin/env bash


#Checking for dependencies

if ! command -v hmmscan > /dev/null; then
        printf "\n HMMER3 not found in path. Exiting"
        exit
fi

if ! command -v efetch > /dev/null; then
	printf "\n Entrez-Direct tools not found in path. Exiting"
	exit
fi

if ! command -v faidx > /dev/null; then
        printf "\n Pyfaidx not found in path. Exiting"
        exit
fi

if ! command -v taxonkit > /dev/null; then
	printf "\n Taxonkit not found in path." 
	printf "\n Activating the included executable"
	chmod +x taxonkit 
fi

if ! command -v FastTree > /dev/null; then
	printf "\n Fasttree not found in path." 
	printf "\n Activating the included executable"
	chmod +x FastTreeDbl 
fi

if ! command -v muscle > /dev/null; then
	printf "\n Muscle not found in path." 
	printf "\n Activating the included executable"
	chmod +x muscle
fi


#Binomica logo

echo "                                                  "
echo "                                                  "
echo "                                                  "
echo "            ,@##          %..      %#/            "
echo "             &%/          %,.      %#/            "
echo "            *&%.          %..      %#(            "
echo "            (%%           %..      %#(            "
echo "            /%%           %*.      %#(            "
echo "             %%            &#.,*,(###.            "
echo "             #%                %,(                "
echo "             #%/               %,(                "
echo "            .%%(               %,/                "
echo "             #%%%%#%%%#%(./(*,**#(                "
echo "                       #*,                        "
echo "                       .%,.                       "
echo "                        %,#                       "
echo "                        %,%                       "
echo "                        %#%                       "
echo "                       (%%%                       "
echo "                        (%                        "
echo "                                                  "
echo "         Binomica Labs - TelePhage v. 0.01        "
echo "                                                  "
echo "             Small Thoughtful Science             "
echo "                                                  "
echo "                                                  "
echo "                                                  "
echo "##################################################"


#Checking for output directory and creating one

if [ -d "./TP_output" ]
then
        echo "TP_output directory already exists. Exiting for now."
        exit
else
        mkdir TP_output
fi

#Checking for Taxonkit dependency directry .taxonkit/, creating one if not found

if [ -d ~/.taxonkit/ ]
then
        echo "                                               "
        echo "                                               "
        echo "   Taxonkit dump directory found. Processing.  "
        echo "                                               "
        echo "                                               "
else
        echo "                                               "
        echo "                                               "
        echo "   Downloading and moving taxonkit dump files. "
        echo "                                               "
        echo "                                               "

        mkdir ~/.taxonkit/
        wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
        tar -xvf taxdump.tar.gz
        mv delnodes.dmp merged.dmp names.dmp nodes.dmp ~/.taxonkit/
        rm citations.dmp division.dmp gc.prt gencode.dmp readme.txt taxdump.tar.gz
fi

#Checking for Pfam.hmm and other required files for Hmmscan step
#If the files are absent the script will attempt to download and hmmpress the file

FILE=Pfam-A.hmm

if [ -f "$FILE" ]; then
	echo "                                           "
	echo "                                           "
	echo "   $FILE found, proceeding with analysis   "
	echo "                                           "
	echo "                                           "
else
	echo "                                                                           "
        echo "                                                                           "	
	echo "   $FILE not found. Downloading matching version and activating hmmpress   "
	echo "                                                                           "
	echo "                                                                           "
	wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam34.0/Pfam-A.hmm.gz
	gunzip Pfam-A.hmm.gz
	hmmpress Pfam-A.hmm
fi

#Finding Pfam accession for input

input_protein="$1"
hmmscan --tblout "$1".pf Pfam-A.hmm "$1"

#Isolating the top scoring pfam accession number from HMMscan result

cat *.pf | grep -i "PF" | awk '{print $2}' | head -n 1 > "$1"_input_protein_pfam.txt 
cp "$1"_input_protein_pfam.txt *.pf TP_output/
rm *.pf

#Screen preformatted phage pfam output for similar protein acccession numbers

input_pfam=$(<"$1"_input_protein_pfam.txt)

echo "                                                  "
echo "##################################################"
echo "##  Top pfam hit for the input is $input_pfam    ##"
echo "##################################################"
echo "                                                  "

zcat allphage-pf34.pfam.gz | grep -i $input_pfam | awk '{print $2}' > "$1"_similar_protein_accessions.acc && echo "Protein accessions list created" || echo "Protein accession list compilation failed" 
cp "$1"_similar_protein_accessions.acc TP_output/

#Extracting full protein sequences associated with the accession numbers
#We're first parsing the accession groups into 200
#Multiple simultaneous download from NCBI is not recommended!

JoinIntoGroupsOf() {
	xargs -n "$@" echo |
	sed 's/ /,/g'
}

cat "$1"_similar_protein_accessions.acc | JoinIntoGroupsOf 200 | xargs -n 1 sh -c \
	'efetch -email "sung@binomicalabs.org" -db protein -id $0 -format fasta' \
	> "$1"_similar_proteins.seq && echo "Protein FASTA for matching accessions compiled" || \
	echo "Protein FASTA compilation failed"
	cp "$1"_similar_proteins.seq TP_output/

cat "$1"_similar_protein_accessions.acc | JoinIntoGroupsOf 200 | xargs -n 1 sh -c \
	'efetch -email "sung@binomicalabs.org" -db protein -id $0 -format fasta_cds_na' \
	> "$1"_similar_proteins_cds.seq && echo "Nucleotide FASTA for matching accessions compiled" || \
	echo "Nucleotide FASTA compilation failed"
	cp "$1"_similar_proteins_cds.seq TP_output/

epost -input "$1"_similar_protein_accessions.acc -db protein | efetch -format gb -mode xml | xtract -pattern \
	GBSeq -element GBSeq_accession-version -block GBQualifier -if GBQualifier_name \
	-equals product -element GBQualifier_value > \
	"$1"_similar_proteins_annotation.txt && echo "Protein accessions annotation file compiled" || \
	echo "Protein accessions annotation compilation failed"
	cp "$1"_similar_proteins_annotation.txt TP_output/

#Quick warning message for the user before the more time consuming processes begin
	
echo "                                                                   "
echo "                                                                   "
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo "~~( Stay awhile and listen; following steps will take some time )~~"
echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo "                                                                   "
echo "                                                                   "


echo "                                                  "
echo "##################################################"
echo "## Aligning the pfam matching protein sequences ##"
echo "##################################################"
echo "                                                  "

#Using muscle to align the pfam sharing protein sequences
#Choosing to use included binary or system call

if ! command -v muscle > /dev/null
then
	./muscle -in "$1"_similar_proteins.seq > "$1"_similar_proteins.aln
else
	muscle -in "$1"_similar_proteins.seq > "$1"_similar_proteins.aln
fi

cp "$1"_similar_proteins.aln TP_output/

echo "                                                  "
echo "##################################################"
echo "##   Running FastTree on the protein alignment  ##"
echo "##################################################"
echo "                                                  "

#Using FastTreeDbl to generate a default phylogenetic tree from the alignment
#Choosing to use included binary or system call

if ! command -v FastTree > /dev/null
then
	./FastTreeDbl  "$1"_similar_proteins.aln > "$1"_similar_proteins.tre
else
	FastTree "$1"_similar_proteins.aln > "$1"_similar_proteins.tre
fi

cp "$1"_similar_proteins.tre TP_output/ 

#Screen output for the taxon id search delay"

echo "                                                  "
echo "##################################################"
echo "##             Crunching taxon data             ##"
echo "##################################################"
echo "                                                  "

#Linking the accession numbers to NCBI uid

zgrep -f "$1"_similar_protein_accessions.acc phage-prot.accession2taxid.gz | awk '{print $2"\t"$3}' > "$1"_similar_proteins_to_uids && echo "Protein accessions to taxonomic ID list completed" || echo "Protein accessions to taxonomic ID list failed"
cp "$1"_similar_proteins_to_uids TP_output/

#Screening input_proteins_to_uids file to isolate uid for taxonomy tagging

cat "$1"_similar_proteins_to_uids | awk '{print $2}' > "$1"_phages_with_similar_proteins.uid && echo "Phage taxonomic ID list completed" || echo "Phage taxonomic ID list compilation failed"
cp "$1"_phages_with_similar_proteins.uid TP_output/

#Using Taxonkit to link NCBI uid with taxonomic names for identified phages
#Also creating a temporary linear uid file to create taxonomic tree list

taxonkit lineage "$1"_phages_with_similar_proteins.uid > "$1"_phages-taxonomic-full.txt && echo "Phage taxonomic name list completed" || echo "Phage taxonomic name compilation failed"
taxonkit lineage -t "$1"_phages_with_similar_proteins.uid > "$1"_phages-taxonomic-lineage.txt && echo "Phage taxonomic lineage list completed" || echo "Phage taxonomic lineage compilation failed"

cp "$1"_phages-taxonomic-full.txt "$1"_phages-taxonomic-lineage.txt TP_output/


#Clean up for the working directory

rm "$1"_input_protein_pfam.txt "$1"_similar_protein_accessions.acc "$1"_similar_proteins.seq "$1"_similar_proteins_cds.seq "$1"_similar_proteins_annotation.txt "$1"_similar_proteins_to_uids "$1"_phages_with_similar_proteins.uid "$1"_phages-taxonomic-full.txt "$1"_phages-taxonomic-lineage.txt "$1"_similar_proteins.aln "$1"_similar_proteins.tre && echo "All TelePhage processes completed. Have a nice day!" || echo "TelePhage process terminated with error" 

echo "                                                  "
echo "                                                  "
echo "                                                  "
