import protein_complex_maps.protein_complex_maps.complex_map_website.complex_db as cdb
#from flask import Flask
from sqlalchemy.orm.exc import NoResultFound
from sqlalchemy import func, or_

import pandas as pd
from flask import Flask, Response, send_file, session

import mpmath as mpm
#import scipy.misc as misc
import StringIO as sio
import tempfile as tf
import datetime as dt

db = cdb.get_db()
app = cdb.get_app()
db.create_all()

#from flask.ext.wtf import Form
from flask_wtf import Form
from wtforms.fields import StringField, SubmitField, TextAreaField

#kdrew: from http://stackoverflow.com/questions/33468821/is-scipy-misc-comb-faster-than-an-ad-hoc-binomial-computation
def choose(n, k):
    """
    A fast way to calculate binomial coefficients by Andrew Dalke (contrib).
    """
    if 0 <= k <= n:
        ntok = 1
        ktok = 1
        for t in xrange(1, min(k, n - k) + 1):
            ntok *= n
            ktok *= t
            n -= 1
        return ntok // ktok
    else:
        return 0

#kdrew: hypergeometric test, k = overlap, n = total number of genes input, m = total number of genes in complex, N = total number of genes in all complexes
def pval(k,n,m,N):
    pv = 0.0                                                 
    #N_choose_m = 1.0*misc.comb(N,m)
    for i in range(k,int(min(m,n)+1)):
        pi = ( mpm.binomial(n,i) * mpm.binomial((N-n), (m-i)) ) / mpm.binomial(N,m)
        #pi = ( misc.comb(n,i) * misc.comb((N-n), (m-i)) ) /  N_choose_m
        #kdrew: somethings wrong with this implementation, returns 0 always (or the cases I've tried)
        #pi = ( choose(n,i) * choose((N-n), (m-i)) ) / choose(N,m)
        pv += pi
    return pv

#kdrew: convert results into csv format
def complexes_to_csv(complexes, prot_ids=[], pvalue_dict=None, genename_cannotfind_errors=None, genename_nocomplex_errors=None):
    complexes_dict = dict()
    complexes_dict['humap_id'] = []
    complexes_dict['rank'] = []
    if pvalue_dict != None:
        complexes_dict['pvalue'] = []
    if len(prot_ids) > 0:
        complexes_dict['searched_accs'] = []
        complexes_dict['searched_genenames'] = []
    complexes_dict['accs'] = []
    complexes_dict['genenames'] = []

    for comp in complexes:
        complexes_dict['humap_id'].append(comp.humap_id)
        complexes_dict['rank'].append(comp.top_rank)
        if pvalue_dict != None:
            complexes_dict['pvalue'].append(pvalue_dict[comp])
        if len(prot_ids) > 0:
            complexes_dict['searched_accs'].append(' '.join([prot.uniprot_acc for prot in comp.proteins if prot.id in prot_ids]))
            complexes_dict['searched_genenames'].append(' '.join([prot.genename() for prot in comp.proteins if prot.id in prot_ids]))
        complexes_dict['accs'].append(' '.join([prot.uniprot_acc for prot in comp.proteins]))
        complexes_dict['genenames'].append(' '.join([prot.genename() for prot in comp.proteins]))
            

    if genename_cannotfind_errors != None:
        complexes_dict['humap_id'].append("COULD_NOT_FIND_GENENAMES")
        complexes_dict['rank'].append(None)
        if pvalue_dict != None:
            complexes_dict['pvalue'].append(None)
        if len(prot_ids) > 0:
            complexes_dict['searched_accs'].append(None)
            complexes_dict['searched_genenames'].append(' '.join(genename_cannotfind_errors))
        complexes_dict['accs'].append(None)
        complexes_dict['genenames'].append(None)

    if genename_nocomplex_errors != None:
        complexes_dict['humap_id'].append("NO_COMPLEXES_FOR_GENENAMES")
        complexes_dict['rank'].append(None)
        if pvalue_dict != None:
            complexes_dict['pvalue'].append(None)
        if len(prot_ids) > 0:
            complexes_dict['searched_accs'].append(None)
            complexes_dict['searched_genenames'].append(' '.join(genename_nocomplex_errors))
        complexes_dict['accs'].append(None)
        complexes_dict['genenames'].append(None)

    df = pd.DataFrame.from_dict(complexes_dict)
    return df.to_csv()


class SearchForm(Form):
    complex_id = StringField(u'Complex ID:')
    #genename = StringField(u'Gene Name (ex. OFD1):')
    listOfGenenames = TextAreaField(u'List of Gene Names or Uniprot ACCs (ex. OFD1 PCM1 Q1MSJ5):')
    enrichment = StringField(u'Enrichment (ex. cilium):')
    protein = StringField(u'Protein (ex. Centrosomal protein):')
    submit = SubmitField(u'Search')

from flask import render_template
from flask import url_for, redirect, request, jsonify

@app.route("/")
def root(complexes=[]):
    return render_template('frontpage.html')

@app.route("/h3")
def show_complexes(complexes=[]):
    print complexes
    #complexes = cdb.Complex.query.all()
    form = SearchForm()
    return render_template('show_complexes.html', form=form, complexes=complexes)
    #return render_template('show_complexes.html', form=form, complexes=complexes, prot_ids=[], pvalue_dict=dict(), error=error)


@app.route("/getComplexesCSV/")
def getComplexesCSV():
    complex_csv_entry = db.session.query(cdb.SessionComplexCSV).filter(cdb.SessionComplexCSV.remote_addr == session['remote_addr'], cdb.SessionComplexCSV.query_id == session['query_id']).first()
    return Response(
            complex_csv_entry.complex_csv,
            mimetype="text/csv",
            headers={"Content-disposition":
                     "attachment; filename=complexes.csv"})

@app.route("/displayComplexesForListOfGeneNames")
def displayComplexesForListOfGeneNames():
    listOfGenenames = request.args.get('listOfGenenames').split()
    form = SearchForm()
    #kdrew: do error checking
    error=None
    genename_cannotfind_errors = []
    genename_nocomplex_errors = []

    #print listOfGenenames

    session['query_id'] = str(dt.datetime.now())
    session['remote_addr'] = request.remote_addr


    all_genes = []
    complexes = []
    all_proteins = []
    error_proteins = []
    error_uniprot_proteins = []

    for genename in listOfGenenames:
        #print genename
        #kdrew: tests to see if genename is a valid genename
        genes = db.session.query(cdb.Gene).filter((func.upper(cdb.Gene.genename) == func.upper(genename))).all()

        uniprot_acc_proteins = db.session.query(cdb.Protein).filter((func.upper(cdb.Protein.uniprot_acc) == func.upper(genename))).all()
        all_proteins = all_proteins + uniprot_acc_proteins

        if len(genes) == 0 and len(uniprot_acc_proteins) == 0:
            #kdrew: input genename is not valid, flash message
            genename_cannotfind_errors.append(genename)

        #all_genes = all_genes + genes

        #print [g.genename for g in all_genes]

        found_complexes_flag = False
        for gene in genes:
            print "current genename: %s" % gene.genename
            try:
                #kdrew: why am I matching based on gene_id and not gene.protein_key == Protein.id, switched and can't see a problem
                #proteins = db.session.query(cdb.Protein).filter((cdb.Protein.gene_id == gene.gene_id)).all()
                proteins = db.session.query(cdb.Protein).filter((cdb.Protein.id == gene.protein_key)).all()

            except NoResultFound:
                if error == None:
                    error = ""
                #kdrew: input genename is not valid, flash message
                genename_cannotfind_errors.append(gene.genename)

            all_proteins = all_proteins + proteins

            error_proteins_current = []
            for protein in proteins:
                print "current protein: %s" % protein.genename()
                print "current protein complexes: %s" % len(protein.complexes)
                print "found_complexes_flag: %s" % found_complexes_flag
                if len(protein.complexes) == 0 and not found_complexes_flag:
                    error_proteins_current.append(protein)
                else:
                    try:
                        complexes = complexes + protein.complexes
                        #kdrew: if found complexes for any of the proteins attached to gene then do not report error
                        error_proteins_current = []
                        found_complexes_flag = True
                    except NoResultFound:
                        continue

            #kdrew: add proteins from this iteration that did not have complexes to the list of error proteins
            error_proteins = error_proteins + error_proteins_current
            print "error_proteins: %s" % ' '.join([p.genename() for p in error_proteins])

        #kdrew: making a small assumption here that genenames and uniprot accs do not overlap, I can't image that they would ever
        found_complexes_flag = False
        error_uniprot_proteins_current = []
        for uniprot_acc_protein in uniprot_acc_proteins:
            if len(uniprot_acc_protein.complexes) == 0 and not found_complexes_flag:
                error_uniprot_proteins_current.append(uniprot_acc_protein)
            else:
                try:
                    complexes = complexes + uniprot_acc_protein.complexes
                    #kdrew: if found complexes for any of the proteins attached to gene then do not report error
                    error_uniprot_proteins_current = []
                    found_complexes_flag = True
                except NoResultFound:
                    continue

        #kdrew: add proteins from this iteration that did not have complexes to the list of error proteins
        error_uniprot_proteins = error_uniprot_proteins + error_uniprot_proteins_current

    #if len(complexes) == 0:
    #    if error == None:
    #        error = ""
    #    error = error + "No complexes found for genenames: %s<br>" % ', '.join([g.genename for g in all_genes])
    #    genename_nocomplex_errors = genename_nocomplex_errors + [g.genename for g in all_genes]

    if len(error_proteins) > 0:
        if error == None:
            error = ""
        error = error + "No complexes found for genenames: %s<br>" % ', '.join([p.genename() for p in error_proteins])
        genename_nocomplex_errors = genename_nocomplex_errors + [p.genename() for p in error_proteins]

    if len(error_uniprot_proteins) > 0:
        if error == None:
            error = ""
        error = error + "No complexes found for genenames: %s<br>" % ', '.join([p.uniprot_acc for p in error_uniprot_proteins])
        genename_nocomplex_errors = genename_nocomplex_errors + [p.uniprot_acc for p in error_uniprot_proteins]


    n = len(all_proteins)
    N = db.session.query(cdb.ProteinComplexMapping).distinct(cdb.ProteinComplexMapping.protein_key).group_by(cdb.ProteinComplexMapping.protein_key).count()
    pvalue_dict = dict()
    pvalue_forsorting_dict = dict()
    for c in set(complexes):
        k = complexes.count(c)
        m = len(c.proteins) 
        print "complex: %s k: %s n: %s m: %s N: %s" % (c.complex_id,k,n,m,N)
        pvalue = pval(k=k,n=n,m=m,N=N)
        pvalue_dict[c] = '{:.2e}'.format(float(pvalue))
        pvalue_forsorting_dict[c] = pvalue
        

    #complexes = list(set(complexes))
    #complexes = [x[1] for x in sorted(((complexes.count(e), e) for e in set(complexes)), reverse=True)]
    #complexes = [x[1] for x in sorted((((complexes.count(e), -1*e.top_rank), e) for e in set(complexes)), reverse=True)]
    complexes = [x[1] for x in sorted((((pvalue_forsorting_dict[e], e.top_rank), e) for e in set(complexes)), reverse=False)]


    
    complex_csv = complexes_to_csv(complexes=complexes, prot_ids=[p.id for p in all_proteins], pvalue_dict=pvalue_dict, genename_cannotfind_errors=genename_cannotfind_errors, genename_nocomplex_errors=genename_nocomplex_errors)
    #kdrew: store in database request.remote_addr
    cdb.SessionComplexCSV().delete_expired_entries()
    complex_csv_entry = cdb.get_or_create(db, cdb.SessionComplexCSV, 
                                remote_addr = session['remote_addr'],
                                query_id = session['query_id'],
                                complex_csv = complex_csv
                                )
    db.session.add(complex_csv_entry)
    db.session.commit()

    #print [p.id for p in all_proteins]
    return render_template('show_complexes.html', form=form, complexes=complexes, prot_ids=[p.id for p in all_proteins], pvalue_dict=pvalue_dict, error=error, genename_cannotfind_errors=genename_cannotfind_errors, genename_nocomplex_errors=genename_nocomplex_errors) 

@app.route("/displayComplexesForEnrichment")
def displayComplexesForEnrichment():
    enrichment = request.args.get('enrichment')
    form = SearchForm()
    error=None

    session['query_id'] = str(dt.datetime.now())
    session['remote_addr'] = request.remote_addr

    #print enrichment
    #kdrew: do error checking
    try:
        enrichment_complex_keys_query = db.session.query(cdb.ComplexEnrichment.complex_key).filter(
                    ( cdb.ComplexEnrichment.description.like('%'+enrichment+'%')) | (cdb.ComplexEnrichment.name.like('%'+enrichment+'%' ) ) )
        enrichment_complex_keys = enrichment_complex_keys_query.all()
        enrichment_complex_keys_set = set([x[0] for x in enrichment_complex_keys])
        if len(enrichment_complex_keys_set) == 0:
            complexes = []
        else:
            complexes = db.session.query(cdb.Complex).filter(cdb.Complex.id.in_(enrichment_complex_keys_set)).all()
    except NoResultFound:
        complexes = []

    if len(complexes) == 0:
        error = "No complexes found for given enrichment term: %s" % enrichment

    complexes = [x[1] for x in sorted((((complexes.count(e), -1*e.top_rank), e) for e in set(complexes)), reverse=True)]

    complex_csv = complexes_to_csv(complexes=complexes, prot_ids=[], pvalue_dict=None, genename_cannotfind_errors=None, genename_nocomplex_errors=None)
    #kdrew: store in database request.remote_addr
    cdb.SessionComplexCSV().delete_expired_entries()
    complex_csv_entry = cdb.get_or_create(db, cdb.SessionComplexCSV, 
                                remote_addr = session['remote_addr'],
                                query_id = session['query_id'],
                                complex_csv = complex_csv
                                )
    db.session.add(complex_csv_entry)
    db.session.commit()

    #return render_template('show_complexes.html', form=form, complexes=complexes, error=error)
    return render_template('show_complexes.html', form=form, complexes=complexes, prot_ids=[], pvalue_dict=None, error=error)

@app.route("/displayComplexesForProtein")
def displayComplexesForProtein():
    protein_search = request.args.get('protein')
    form = SearchForm()
    error = None

    session['query_id'] = str(dt.datetime.now())
    session['remote_addr'] = request.remote_addr

    #print protein
    #kdrew: do error checking
    complexes = []
    try:
        proteins = db.session.query(cdb.Protein).filter(cdb.Protein.proteinname.like('%'+protein_search+'%')).all()
        for protein in proteins:
            print protein
            complexes = complexes + protein.complexes

        #kdrew: remove redudant complexes
        complexes = list(set(complexes))

    except NoResultFound:
        complexes = []

    if len(complexes) == 0:
        error = "No complexes found for given search term: %s" % protein_search

    complexes = [x[1] for x in sorted((((complexes.count(e), -1*e.top_rank), e) for e in set(complexes)), reverse=True)]

    complex_csv = complexes_to_csv(complexes=complexes, prot_ids=[], pvalue_dict=None, genename_cannotfind_errors=None, genename_nocomplex_errors=None)
    #kdrew: store in database request.remote_addr
    cdb.SessionComplexCSV().delete_expired_entries()
    complex_csv_entry = cdb.get_or_create(db, cdb.SessionComplexCSV, 
                                remote_addr = session['remote_addr'],
                                query_id = session['query_id'],
                                complex_csv = complex_csv
                                )
    db.session.add(complex_csv_entry)
    db.session.commit()

    #return render_template('show_complexes.html', form=form, complexes=complexes, error=error)
    return render_template('show_complexes.html', form=form, complexes=complexes, prot_ids=[], pvalue_dict=None, error=error)

@app.route("/displayComplexes")
def displayComplexes():
    complex_key = request.args.get('complex_key')
    form = SearchForm()
    error=None
    #kdrew: do error checking
    try:
        #comp = db.session.query(cdb.Complex).filter_by(complex_id=complex_key).one()
        comp = db.session.query(cdb.Complex).filter_by(humap_id=complex_key).one()
    except NoResultFound:
        comp = None

    if comp == None:
        error = "No complexes found: %s" % complex_key

    return render_template('complex.html', form=form, comp=comp, error=error)


@app.route(u'/search', methods=[u'POST'])
def searchComplexes():
    form = SearchForm()
    complexes = []
    if form.validate_on_submit():
        #if len(form.genename.data) > 0:
        #    return redirect(url_for('displayComplexesForGeneName', genename=form.genename.data))
        #elif len(form.listOfGenenames.data) > 0:
        if len(form.listOfGenenames.data) > 0:
            return redirect(url_for('displayComplexesForListOfGeneNames', listOfGenenames=form.listOfGenenames.data))
        elif len(form.enrichment.data) > 0:
            return redirect(url_for('displayComplexesForEnrichment', enrichment=form.enrichment.data))
        elif len(form.protein.data) > 0:
            return redirect(url_for('displayComplexesForProtein', protein=form.protein.data))


    #kdrew: added hoping it would fix redirect problem on stale connections
    #return render_template('show_complexes.html', form=form, complexes=complexes)
    return render_template('show_complexes.html', form=form, complexes=complexes, prot_ids=[], pvalue_dict=None, error="")

@app.route("/about")
def displayAbout():
    return render_template('about.html')

@app.route("/download")
def displayDownload():
    return render_template('download.html')

#kdrew: frontpage.html is the main page for proteincomplexes.org that links to pages of other maps (i.e. humap2, rnaMAP, etc)
@app.route("/frontpage")
def displayFrontPage():
    return render_template('frontpage.html')

#kdrew: about_proteincomplexes.html is the main page for proteincomplexes.org that links to pages of other maps (i.e. humap2, rnaMAP, etc)
@app.route("/about_proteincomplexes")
def displayAboutProteinComplexes():
    return render_template('about_proteincomplexes.html')


if __name__ == "__main__":
    db.create_all()  # make our sqlalchemy tables
    app.run(threaded=True, host='0.0.0.0', port=5000)

#@app.route('/')
#def hello_world():
#   return "Hello, World!"
#
#if __name__ == "__main__":
#   app.run(host='0.0.0.0', port=5000)


