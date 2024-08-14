
import argparse
import numpy as np
import itertools as it
import pandas as pd

import csv

import protein_complex_maps.complex_map_website.complex_db as cdb

def main():

    parser = argparse.ArgumentParser(description="Outputs complex database edges to flat files")
    parser.add_argument("--output_file", action="store", dest="output_file", required=True, 
                                    help="Filename of output file")
    parser.add_argument("--cytoscape", action="store_true", dest="cytoscape", required=False, default=False, 
                                    help="Format output to be read in by cytoscape")

    args = parser.parse_args()

    db = cdb.get_db()
    app = cdb.get_app()

    db.create_all()

    fout = open(args.output_file,"w")
    if args.cytoscape:
        fout.write("shared_name,score,protHD2_score\n")
    else:
        fout.write("uniprot_acc1,uniprot_acc2,score,protHD2_score\n")
    for c in db.session.query(cdb.Complex).all():
        if c:
            if len(c.edges) == 0:
                continue
            print("complex id: %s" % c.complex_id)
            cytoscape_format_proteins = ""
            for e in c.edges:
                print("score: %s" % e.score)
                print("proteins: %s" % ' '.join([p.uniprot_acc for p in e.get_proteins()]))
                p1 = None
                p2 = None
                for p in e.get_proteins():
                    if p1 == None:
                        p1 = p.uniprot_acc
                    else:
                        p2 = p.uniprot_acc
                cytoscape_format_proteins = ' ((pp)) '.join(["%s_%s" % (c.complex_id, p) for p in [p1,p2]])
                #kdrew: reverse protein names so both entries will be available for matching in cytoscape
                cytoscape_format_proteins_rev = ' ((pp)) '.join(["%s_%s" % (c.complex_id, p) for p in [p2,p1]])
                print(cytoscape_format_proteins)
                if args.cytoscape:
                    fout.write("%s,%s,%s\n" % (cytoscape_format_proteins, e.score, e.get_prothd_score()))
                    fout.write("%s,%s,%s\n" % (cytoscape_format_proteins_rev, e.score, e.get_prothd_score()))
                else:
                    fout.write("%s,%s,%s\n" % (','.join([p.uniprot_acc for p in e.get_proteins()]), e.score, e.get_prothd_score()))
        else:
            print("Cannot find complex %s" % (c.complex_id))

    fout.close()

if __name__ == "__main__":
    main()


