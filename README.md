# Telephage 

### From proteins of interest to phage phylogeny in single line!


            ,@##          %..      %#/            
             &%/          %,.      %#/            
            *&%.          %..      %#(            
            (%%           %..      %#(            
            /%%           %*.      %#(            
             %%            &#.,*,(###.            
             #%                %,(                
             #%/               %,(                
            .%%(               %,/                
             #%%%%#%%%#%(./(*,**#(                
                       #*,                        
                       .%,.                       
                        %,#                       
                        %,%                       
                        %#%                       
                       (%%%                       
                        (%                                                                          "
                 
           **Binomica Labs - TelePhage**         
              
            *Small Thoughtful Science*             

---

Telephage is a simple tool meant to immediately provide a library of interesting data linking submitted proteins to known phages. 
The goal of the package is to help beginning amateur researchers gain a broader perspective on any particular protein by focusing on immediate screening of known databases with an eye for possible evolutionary relationships.

Simply clone or download the folder, and execute the Telephage.sh script followed by a protein fasta file of your choice. 

`./Telephage.sh myown.faa`

The script will download and setup the necessary Pfam and Taxonkit dump file if necessary, and point out any missing dependencies. 

Once finished, the script will produce a folder named TP_output containing about 12 files including:

- Pfam screening result of the submitted protein
- Genbank accession list of all ncbi phage proteins that belongs to the same Pfam group
- Annotation of the phage proteins accessions
- .seq file containing all the protein sequences of the phage proteins
- .seq file containing the CDS nucleotide sequences of the phage proteins, along with their coordinates 
- .aln Muscle alignment of all phage proteins that belong to the Pfam group of your protein of interest
- .tre FastTree cladogram of the similarity/distance between the phage protein sequences 
- .txt files containing taxonomic info and lineage on all phages that contain the proteins processed above

---

**Dependencies**

HMMER and Entrez-Direct are absolute requirements. 

It is highly recommended that you have Muscle, FastTree, and Taxonkit installed on your system path already. 
However, if the script fails to detect any of the required three programs it will use the binaries included in the directory. 





