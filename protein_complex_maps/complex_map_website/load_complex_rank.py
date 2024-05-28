
import argparse
import numpy as np
import itertools as it
import pandas as pd

import csv

import protein_complex_maps.complex_map_website.complex_db as cdb

def main():

    parser = argparse.ArgumentParser(description="Loads complex rank from input files")
    parser.add_argument("--toprank_file", action="store", dest="toprank_file", required=True, 
                                    help="Complex top rank filename")

    args = parser.parse_args()

    db = cdb.get_db()
    app = cdb.get_app()

    db.create_all()

    rank_df = pd.read_csv(args.toprank_file)
    rank_dict = dict()
    for _, row in rank_df.iterrows():

        rank_dict[frozenset(row.Uniprot_ACCs.split())] = row.Cluster_Num

    for c in db.session.query(cdb.Complex).all():
        print(c.complex_id)
        print(len(c.proteins))
        print(c.edges)
        print(c.sorted_proteins())
        pset = frozenset([p.uniprot_acc for p in c.proteins])
        print(pset)
        print(rank_dict[pset])
        c.top_rank = rank_dict[pset]
        db.session.commit()

        #c = db.session.query(cdb.Complex).filter_by(complex_id=row.clust_id).first()
        #if c:
        #    print "complex id: %s" % c.id
        #    c.top_rank = row.top_rank
        #    db.session.commit()
        #else:
        #    print "Cannot find complex %s" % (complex_id)


if __name__ == "__main__":
    main()


