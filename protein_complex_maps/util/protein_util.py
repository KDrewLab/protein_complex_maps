
#kdrew: this file is for utilities for querying features about proteins
from __future__ import print_function

#import urllib, urllib2 
import urllib.request, urllib.parse, urllib.error

import re
import requests
import difflib
import Bio.Seq as Seq
import Bio.SeqRecord as SeqRecord

ACC_QUERY_LENGTH = 100
FUZZY_MATCH_THRESHOLD = 0.25

#claude: UniProt retired the legacy www.uniprot.org/uniprot/?format=tab API (2022); use the REST API stream endpoint (returns all results, no pagination)
UNIPROT_STREAM_URL = "https://rest.uniprot.org/uniprotkb/stream"

#claude: map legacy UniProt column names (still used by callers) to REST API field names; unknown keywords are passed through unchanged
UNIPROT_FIELD_MAP = {
    'id': 'accession',
    'length': 'length',
    'organism': 'organism_name',
    'genes': 'gene_names',
    'genes(PREFERRED)': 'gene_primary',
    'protein+names': 'protein_name',
    'protein names': 'protein_name',
    'sequence': 'sequence',
    'reviewed': 'reviewed',
}

#claude: UniProtKB accession format (https://www.uniprot.org/help/accession_numbers); the REST API rejects a whole query (HTTP 400) if any accession is malformed
UNIPROT_ACC_REGEX = re.compile(r'^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$')

##kdrew: queries uniprot for protein sequence length
#def get_length_uniprot( protein_id ):
#   length_dict = dict()
#   report_query = "http://www.uniprot.org/uniprot/?format=tab&query=accession:%s&columns=id,length" % (protein_id,)
#   f = urllib2.urlopen(report_query)
#   for line in f.readlines():
#       if protein_id == line.split()[0]:
#           return int(line.split()[1])
#   return None


def get_organism_uniprot( protein_ids ):
    return get_from_uniprot( protein_ids, "organism" )

#kdrew: queries uniprot for protein sequence length
def get_length_uniprot( protein_ids ):
    return get_from_uniprot( protein_ids, "length" )

#kdrew: queries uniprot for protein sequence genename
def get_genenames_uniprot( protein_ids ):
    return get_from_uniprot( protein_ids, "genes(PREFERRED)" )

def get_sequences_uniprot( protein_ids, seqrecord=False ):
    seq_map = get_from_uniprot(protein_ids, keyword="sequence") 
    if seqrecord:
        return_seq_map = dict()
        organism_map = get_organism_uniprot( protein_ids )
        for key in seq_map:
            if seq_map[key] != None:
                #claude: Bio.Alphabet was removed in Biopython 1.78; molecule type is now a SeqRecord annotation
                s = Seq.Seq(seq_map[key])
                sr = SeqRecord.SeqRecord(s, id=key, description=organism_map[key], annotations={"molecule_type": "protein"})
                #sr = SeqRecord.SeqRecord(Seq.Seq(seq_map[key], ba.generic_protein), id=key, description=organism_map[key])
                #print(sr)
                return_seq_map[key] = sr
        return return_seq_map

    return seq_map

def get_from_uniprot( protein_ids, keyword, return_list=False ):
    return_dict = dict()

    #kdrew: sometimes uniprot accs have added '-1' that does not play well with this webservice
    cleaned_protein_ids = []
    for prot in protein_ids:
        cleaned_protein_ids.append(prot.split('-')[0])

    protein_ids = protein_ids + cleaned_protein_ids

    #claude: only query unique, well-formed accessions; results are keyed by the base accession, which is in protein_ids as before
    query_ids = []
    for prot in cleaned_protein_ids:
        if prot in query_ids:
            continue
        if UNIPROT_ACC_REGEX.match(prot):
            query_ids.append(prot)
        else:
            print("WARNING: skipping malformed UniProt accession: %s" % prot)

    field = UNIPROT_FIELD_MAP.get(keyword, keyword)

    #kdrew (modified by claude): break up query into chunks so as not to get 414 error
    for i in range( (len(query_ids)//ACC_QUERY_LENGTH)+1 ):
        start_splice = i*ACC_QUERY_LENGTH
        stop_splice = (i+1)*ACC_QUERY_LENGTH
        chunk = query_ids[start_splice:stop_splice]
        if len(chunk) == 0:
            continue
        params = {
            'query': " OR ".join(["accession:%s" % prot for prot in chunk]),
            'fields': "accession,%s" % field,
            'format': 'tsv',
        }
        try:
            f = requests.get(UNIPROT_STREAM_URL, params=params)
            f.raise_for_status()
        except requests.exceptions.HTTPError as msg:
            print(msg)
            continue

        for line in f.text.split('\n'):
            #claude: skip blank lines (response text ends with a newline)
            if len(line.split()) == 0:
                continue
            if line.split()[0] in protein_ids:
                if field == "length":
                    return_dict[line.split()[0]] = int(line.split()[1])
                elif field == "organism_name":
                    return_dict[line.split()[0]] = line.split('\t')[1]
                elif field == "protein_name":
                    return_dict[line.split()[0]] = line.split('\t')[1]
                else:
                    if return_list:
                        try:
                            return_dict[line.split()[0]] = line.split()[1:]
                        except IndexError:
                            return_dict[line.split()[0]] = []
                    else:
                        try:
                            return_dict[line.split()[0]] = line.split()[1]
                        except IndexError:
                            return_dict[line.split()[0]] = None


    return return_dict

def get_from_uniprot_by_genename( gene_ids, organism="", gene_prefix="gene_exact", reviewed=False ):
    return_dict = dict()

    gene_ids_case_incensitive = gene_ids + [x.lower() for x in gene_ids] + [y.upper() for y in gene_ids] 

    reviewed_str = ""
    if reviewed:
        reviewed_str = "+and+reviewed:yes"

    #kdrew: break up query into chunks so as not to get 414 error
    for i in range( (len(gene_ids)//ACC_QUERY_LENGTH)+1 ):
        start_splice = i*ACC_QUERY_LENGTH
        stop_splice = (i+1)*ACC_QUERY_LENGTH                                                                                 
        genes_formatted = ["%s:%s" % (gene_prefix, id1,) for id1 in gene_ids[start_splice:stop_splice]]
        report_query = "https://www.uniprot.org/uniprot/?format=tab&query=organism:(%s)+and+(%s)%s" % (organism, "+or+".join(genes_formatted), reviewed_str )
        print(report_query)
        try:
            #f = urllib2.urlopen(encoded_query)
            f = requests.get(report_query)
            #print(f.text)
            for line in f.text.split('\n'):
                print(line)
                #kdrew: the returned line has a bunch of extra gene names (and other info), finding the intersection of the gene ids and the line returns the original queried gene id
                intersection = list(set(line.split()).intersection(gene_ids_case_incensitive))
                if intersection:
                    try:
                        return_dict[intersection[0]].append(line.split()[0])
                    except:
                        return_dict[intersection[0]] = [line.split()[0],] 
        except urllib.error.HTTPError as e:
            print(e) #, e.read()
            continue


    return return_dict


#kdrew: uses uniprot webservice to map ids
#kdrew: from_id and to_id are abbreviations of dbid names which can be found: http://www.uniprot.org/faq/28
def map_protein_ids( id_list, from_id, to_id, contact_email, reviewed=False ):
    url = 'https://www.uniprot.org/mapping/'

    query_str = ' '.join(id_list)
    params = {
        'query':query_str,
        'from':'%s' % (from_id,),
        'to':'%s' % (to_id,),
        'format':'tab',
        #kdrew: noticed a TypeError exception in join but not sure why
    }

    print(query_str)

    data = urllib.parse.urlencode(params)
    print(data)
    request = urllib.request.Request(url, data)
    contact = contact_email
    request.add_header('User-Agent', 'Python %s' % contact)
    response = urllib.request.urlopen(request)

    #kdrew: put resulting map into dictionary of lists (one to many)
    return_dict = {}
    for line in response.readlines():
        print(line)
        if line.split()[0] != "From":
            try:
                if len(line.split()) > 2:
                    return_dict[line.split()[0]] = return_dict[line.split()[0]] + line.split()[1:]
                else:
                    return_dict[line.split()[0]].append(line.split()[1])
            except:
                if len(line.split()) > 2:
                    return_dict[line.split()[0]] = line.split()[1:]
                else:
                    return_dict[line.split()[0]] = [line.split()[1],] 


    for i in id_list:
        if i not in list(return_dict.keys()):
            print("No id match for %s, adding empty list" % i)
            return_dict[i] = []

    if reviewed:
        #kdrew: flatten ids into set
        ret_id_list = list(set([item for sublist in list(return_dict.values()) for item in sublist]))
        #print ret_id_list
        is_reviewed_list = get_from_uniprot(ret_id_list, 'reviewed')
        #print is_reviewed_list
        unreviewed_list = [item for item in is_reviewed_list if is_reviewed_list[item] == 'unreviewed']
        #print unreviewed_list

        for i in return_dict:
            return_dict[i] = list(set(return_dict[i]) - set(unreviewed_list))

    return return_dict

