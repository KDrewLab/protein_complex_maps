
import argparse
import numpy as np
import pandas as pd
import itertools as it
from sqlalchemy import func, or_, and_

import csv

import protein_complex_maps.complex_map_website.complex_db as cdb

def main():

    parser = argparse.ArgumentParser(description="Loads ProtHD table")
    parser.add_argument("--prothd_file", action="store", dest="prothd_file", required=True, 
                                    help="Filename ProtHD table")

    args = parser.parse_args()

    db = cdb.get_db()
    app = cdb.get_app()

    db.create_all()

    prothd_df = pd.read_csv(args.prothd_file)
    #kdrew: iterate through complexes
    for _, row in prothd_df.iterrows():

        prot1 = db.session.query(cdb.Protein).filter_by(uniprot_acc=row.Protein_1).first()
        prot2 = db.session.query(cdb.Protein).filter_by(uniprot_acc=row.Protein_2).first()
        if prot1 == None or prot2 == None:
            continue
        if prot2.id < prot1.id:
            prot2, prot1 = prot1, prot2

            e = db.session.query(cdb.Edge).filter( and_(cdb.Edge.protein_key == prot1.id, cdb.Edge.protein_key2 == prot2.id) ).first()
            if e != None:
                print("prot1: %s, prot2: %s " % (prot1.id, prot2.id))
                print("edge id: %s" % e.id)
                print("score: %s" % row.RF_covariation_prob)
                ecm = cdb.get_or_create(db, cdb.ProtHD, edge_key=e.id, prothd_score=row.RF_covariation_prob)
                db.session.add(ecm)
                db.session.commit()


if __name__ == "__main__":
    main()


