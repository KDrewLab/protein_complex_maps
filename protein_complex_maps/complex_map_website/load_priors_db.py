
import argparse
import numpy as np
import itertools as it

import pandas as pd
import csv

import protein_complex_maps.complex_map_website.complex_db as cdb

def main():

    #kdrew: should really generalize the bait prey files and make them optional, also need to reorganize the readin through pandas
    #kdrew: but just trying to get the db loaded

    parser = argparse.ArgumentParser(description="Loads edge sql tables from input files")
    parser.add_argument("--edge_file", action="store", dest="edge_file", required=True, 
                                    help="Filename edge table")

    args = parser.parse_args()

    db = cdb.get_db()
    app = cdb.get_app()

    db.create_all()

    BATCH_SIZE = 10000

    df = pd.read_csv(args.edge_file)
    priors_dict = {frozenset(x[0:2]):x[2] for x in df[['ID1','ID2','humap3_Prediction']].values }

    updated = 0
    for e in db.session.query(cdb.Edge).all():
        prots = e.get_proteins()
        acc1 = prots[0].uniprot_acc
        acc2 = prots[1].uniprot_acc
        fset = frozenset([str(acc1),str(acc2)])
        try:
            #print(priors_dict[fset])
            e.humap3_score = priors_dict[fset]
            #db.session.commit()
        except KeyError:
            print("no prior for %s" % fset)
            continue

        updated += 1

        # Periodic commit for huge tables
        if updated % BATCH_SIZE == 0:
            db.session.commit()
            print("Committed %s updates" % updated)

    # Final commit
    db.session.commit()


if __name__ == "__main__":
    main()


