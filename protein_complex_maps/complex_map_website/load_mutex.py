
import argparse
import numpy as np
import pandas as pd
import itertools as it
from sqlalchemy import func, or_, and_

import csv

import protein_complex_maps.complex_map_website.complex_db as cdb

def main():

    parser = argparse.ArgumentParser(description="Loads ProtHD table")
    parser.add_argument("--mutex_file", action="store", dest="mutex_file", required=True, 
                                    help="Filename Mutually Exclusive table")

    args = parser.parse_args()

    db = cdb.get_db()
    app = cdb.get_app()

    db.create_all()

    f = open(args.mutex_file)
    #kdrew: iterate through file
    for line in f.readlines():
        tokens = line.split(',')
        acc1 = tokens[0]
        acc2 = tokens[1]
        interface = tokens[2].strip()

        #kdrew: initialize protein instances as None
        prot1, prot2 = None, None

        prot1 = db.session.query(cdb.Protein).filter_by(uniprot_acc=acc1).first()
        prot2 = db.session.query(cdb.Protein).filter_by(uniprot_acc=acc2).first()

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
            
            print("interface: %s" % interface)
            e.interface_mutex = interface
            db.session.commit()


if __name__ == "__main__":
    main()


