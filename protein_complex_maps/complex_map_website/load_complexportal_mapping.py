
import argparse
import sys

import csv

import protein_complex_maps.complex_map_website.complex_db as cdb

def main():

    parser = argparse.ArgumentParser(description="Loads yeast.MAP to ComplexPortal complex mapping")
    parser.add_argument("--complexportal_mapping_file", action="store", dest="mapping_file", required=True,
                                    help="Filename of tab separated mapping file with header (columns used: yeastMAP_ID, CPX_ID, Match_type)")

    args = parser.parse_args()

    db = cdb.get_db()
    app = cdb.get_app()

    db.create_all()

    complex_dict = dict()
    for c in db.session.query(cdb.Complex).all():
        complex_dict[c.complexmap_id] = c

    rows = list(csv.DictReader(open(args.mapping_file), delimiter='\t'))

    #kdrew: check everything before loading anything
    errors = []
    for row in rows:
        if row['yeastMAP_ID'] not in complex_dict:
            errors.append("cannot find complex: %s" % row['yeastMAP_ID'])
        if row['Match_type'] not in ('exact', 'subset'):
            errors.append("unknown match type: %s %s" % (row['yeastMAP_ID'], row['Match_type']))
        if not row['CPX_ID'].startswith('CPX-'):
            errors.append("unexpected ComplexPortal id: %s %s" % (row['yeastMAP_ID'], row['CPX_ID']))
    if len(errors) > 0:
        for e in errors:
            print(e)
        sys.exit(1)

    #kdrew: replace any previous load
    db.session.query(cdb.ComplexPortalMapping).delete()

    for row in rows:
        db.session.add(cdb.ComplexPortalMapping(complex_key=complex_dict[row['yeastMAP_ID']].id, complexportal_id=row['CPX_ID'], match_type=row['Match_type']))
    db.session.commit()

    print("loaded %s mappings for %s complexes" % (len(rows), len(set([row['yeastMAP_ID'] for row in rows]))))
    unmapped = [cid for cid in complex_dict if cid not in set([row['yeastMAP_ID'] for row in rows])]
    print("complexes without a ComplexPortal mapping: %s" % len(unmapped))

if __name__ == "__main__":
    main()
