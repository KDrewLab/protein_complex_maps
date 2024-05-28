
import argparse
import numpy as np
import itertools as it
import pandas as pd

import csv

import protein_complex_maps.complex_map_website.complex_db as cdb

def main():

    parser = argparse.ArgumentParser(description="Loads complex humap_id from input files")
    parser.add_argument("--humap_id_file", action="store", dest="humap_id_file", required=True, 
                                    help="Complex humap_id filename")

    args = parser.parse_args()

    db = cdb.get_db()
    app = cdb.get_app()

    db.create_all()

    id_df = pd.read_csv(args.humap_id_file)
    id_dict = dict()
    for _, row in id_df.iterrows():

        id_dict[frozenset(row.uniprotACCs.split())] = row.clustID

    for c in db.session.query(cdb.Complex).all():
        print(c.complex_id)
        print(len(c.proteins))
        print(c.edges)
        print(c.sorted_proteins())
        pset = frozenset([p.uniprot_acc for p in c.proteins])
        print(pset)
        print(id_dict[pset])
        c.humap_id = id_dict[pset]
        db.session.commit()

if __name__ == "__main__":
    main()


