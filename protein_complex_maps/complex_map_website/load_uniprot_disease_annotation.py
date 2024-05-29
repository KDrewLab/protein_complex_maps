
import argparse
import numpy as np
import itertools as it
import pandas as pd

import csv

import protein_complex_maps.complex_map_website.complex_db as cdb

def main():

    parser = argparse.ArgumentParser(description="Loads disease annotations from input files")
    parser.add_argument("--disease_file", action="store", dest="disease_file", required=True, 
                                    help="Uniprot disease annotation filename")

    args = parser.parse_args()

    db = cdb.get_db()
    app = cdb.get_app()

    db.create_all()

    disease_df = pd.read_table(args.disease_file)
    for _, row in disease_df.iterrows():

        p = db.session.query(cdb.Protein).filter_by(uniprot_acc=row.Entry).first()
        if p:
            print "protein id: %s" % p.id

            if row.Involvement_in_disease != "":
                diseaseAnnotation = cdb.get_or_create(db,cdb.ProteinDiseaseAnnotation, protein_key = p.id, disease_database = 'Uniprot', disease_annotation = row.Involvement_in_disease)
                db.session.commit()

        else:
            print "Cannot find protein %s" % (row.Entry)


if __name__ == "__main__":
    main()


