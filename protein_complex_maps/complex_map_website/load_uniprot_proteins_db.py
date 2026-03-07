
import argparse
import numpy as np
import itertools as it

import csv

import protein_complex_maps.complex_map_website.complex_db as cdb


def main():

    parser = argparse.ArgumentParser(description="Loads uniprot entries into database")
    parser.add_argument("--uniprot_file", action="store", dest="uniprot_file", required=True, 
                            help="Uniprot filename (format: ACC Protein_name    GeneID  primary_genename    genenames)")

    args = parser.parse_args()

    db = cdb.get_db()
    app = cdb.get_app()

    db.create_all()

    f = open(args.uniprot_file,"rb")
    for line in f.readlines():
        lsplit = line.split("\t")
        ACC = lsplit[0]
        entryName = lsplit[2]
        pname = lsplit[3]
        genenames = lsplit[4].split()
        uniprot_url = "http://www.uniprot.org/uniprot/%s" % (ACC)


        #kdrew: some entries in uniprot annotation file do not have geneid so cannot map protein in database exactly, just search based on uniprot acc instead
        #if len(geneids) == 1:
        protein = db.session.query(cdb.Protein).filter_by(uniprot_acc=ACC).first()
        if protein:
            if len(entryName) > 0: 
                gene = cdb.get_or_create(db,cdb.Gene, genename=entryName, protein_key = protein.id)
                db.session.add(gene)
                db.session.commit()

            for genename in genenames:
                #kdrew: do a little clean up
                gname = genename.split(';')[0]
                gene = cdb.get_or_create(db,cdb.Gene, genename=gname, protein_key = protein.id)
                db.session.add(gene)
                db.session.commit()


if __name__ == "__main__":
    main()


