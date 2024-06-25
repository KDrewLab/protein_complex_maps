
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

        #kdrew: initialize protein instances as None
        prot1, prot2 = None, None
        #kdrew: iterate through all protein ids (; separated)
        for acc in row.Protein_1.split(';'):
            print("acc1: %s" % acc)
            prot1 = db.session.query(cdb.Protein).filter_by(uniprot_acc=acc).first()
            print("prot1 id: %s" % prot1)
            if prot1 != None:
                #kdrew: found a protein entry, break out of for loop
                break
        for acc in row.Protein_2.split(';'):
            print("acc2: %s" % acc)
            prot2 = db.session.query(cdb.Protein).filter_by(uniprot_acc=acc).first()
            print("prot2 id: %s" % prot2)
            if prot2 != None:
                #kdrew: found a protein entry, break out of for loop
                break 

        #kdrew: if no proteins are found, continue through loop
        if prot1 == None or prot2 == None:
            continue

        if prot2.id < prot1.id:
            prot2, prot1 = prot1, prot2

        print("Prior to edge look up: prot1: %s, prot2: %s " % (prot1.id, prot2.id))
        e = db.session.query(cdb.Edge).filter( and_(cdb.Edge.protein_key == prot1.id, cdb.Edge.protein_key2 == prot2.id) ).first()
        print e
        if e != None:
            print("prot1: %s, prot2: %s " % (prot1.id, prot2.id))
            print("edge id: %s" % e.id)
            print("score: %s" % row.RF_covariation_prob)
            #kdrew: temp comment
            ecm = cdb.get_or_create(db, cdb.ProtHD, edge_key=e.id, prothd_score=row.RF_covariation_prob)
            db.session.add(ecm)
            db.session.commit()


if __name__ == "__main__":
    main()


